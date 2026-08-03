#!/usr/bin/env python3

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Peak:
    chrom: str
    start: int
    end: int
    last_value: str | None = None

    @property
    def label(self):
        return f"{self.chrom}:{self.start}-{self.end}"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Cluster nearby BED peaks and score clusters by peak count and spacing. "
            "Two consecutive peaks join the same cluster when their gap is <= "
            "--max-gap. Output clusters are ranked by score from highest to lowest."
        )
    )
    parser.add_argument("input", help="Input BED file with at least three columns.")
    parser.add_argument("output", help="Output BED file for clustered regions.")
    parser.add_argument(
        "--max-gap",
        type=int,
        default=10000,
        help="Maximum gap in bp between consecutive peaks in one cluster (default: 10000).",
    )
    parser.add_argument(
        "--score-decay",
        type=float,
        default=None,
        help=(
            "Distance-decay scale in bp for scoring (default: --max-gap). "
            "Score = peak_count * exp(-mean_gap / score_decay)."
        ),
    )
    parser.add_argument(
        "--min-peaks",
        type=int,
        default=3,
        help="Minimum number of peaks required for output (default: 3).",
    )
    parser.add_argument(
        "--stats-tsv",
        default="brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.tsv",
        help=(
            "Peak statistics TSV containing peak_id, log2FoldChange_vs_all_rest, "
            "and padj_vs_all_rest (default: brca_subtype_vs_all_rest_deseq2/"
            "Basal_specific_regions.tsv)."
        ),
    )
    parser.add_argument(
        "--mean-cpm-tsv",
        default=(
            "data_m.Pancancer_peak.specific.mean_replicates."
            "mean_cpm_by_cancer_type.tsv"
        ),
        help="Peak-by-cancer-type mean CPM matrix.",
    )
    parser.add_argument(
        "--var-cpm-tsv",
        default=(
            "data_m.Pancancer_peak.specific.mean_replicates."
            "var_cpm_by_cancer_type.tsv"
        ),
        help="Peak-by-cancer-type CPM variance matrix.",
    )
    return parser.parse_args()


def chromosome_sort_key(chrom):
    match = re.fullmatch(r"chr(\d+|X|Y|M|MT)", chrom, flags=re.IGNORECASE)
    if not match:
        return (1, chrom)
    label = match.group(1).upper()
    if label.isdigit():
        return (0, int(label))
    return (0, {"X": 23, "Y": 24, "M": 25, "MT": 25}[label])


def load_peaks(path):
    peaks = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 columns")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: start and end must be integers"
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: expected 0 <= start < end, got {start}-{end}"
                )
            last_value = fields[-1] if len(fields) > 3 else None
            peaks.append(Peak(fields[0], start, end, last_value))
    return sorted(
        peaks, key=lambda peak: (chromosome_sort_key(peak.chrom), peak.start, peak.end)
    )


def load_peak_stats(path):
    required_columns = {
        "peak_id",
        "log2FoldChange_vs_all_rest",
        "padj_vs_all_rest",
    }
    stats = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = required_columns - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                f"{path}: missing required columns: {sorted(missing_columns)}"
            )
        for line_number, row in enumerate(reader, start=2):
            peak_id = row["peak_id"]
            if peak_id in stats:
                raise ValueError(f"{path}:{line_number}: duplicate peak_id: {peak_id}")
            try:
                stats[peak_id] = {
                    "log2FoldChange": float(row["log2FoldChange_vs_all_rest"]),
                    "FDR": float(row["padj_vs_all_rest"]),
                }
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid numeric value for peak_id {peak_id}"
                ) from exc
    return stats


def load_peak_matrix_means(path, selected_peak_ids):
    means = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "peak_id":
            raise ValueError(f"{path}: first column must be peak_id")
        value_columns = reader.fieldnames[1:]
        if not value_columns:
            raise ValueError(f"{path}: expected at least one value column")
        for line_number, row in enumerate(reader, start=2):
            peak_id = row["peak_id"]
            if peak_id not in selected_peak_ids:
                continue
            if peak_id in means:
                raise ValueError(f"{path}:{line_number}: duplicate peak_id: {peak_id}")
            try:
                values = [float(row[column]) for column in value_columns]
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid numeric value for peak_id {peak_id}"
                ) from exc
            means[peak_id] = sum(values) / len(values)
    missing_peak_ids = sorted(selected_peak_ids - set(means))
    if missing_peak_ids:
        raise ValueError(
            f"{path}: missing {len(missing_peak_ids)} BED peak IDs, including: "
            + ", ".join(missing_peak_ids[:10])
        )
    return means


def cluster_peaks(peaks, max_gap):
    cluster = []
    cluster_end = None
    for peak in peaks:
        if (
            cluster
            and peak.chrom == cluster[-1].chrom
            and peak.start - cluster_end <= max_gap
        ):
            cluster.append(peak)
            cluster_end = max(cluster_end, peak.end)
            continue
        if cluster:
            yield cluster
        cluster = [peak]
        cluster_end = peak.end
    if cluster:
        yield cluster


def summarize_cluster(cluster, score_decay, peak_stats, peak_mean_cpm, peak_var_cpm):
    gaps = [
        max(0, current.start - previous.end)
        for previous, current in zip(cluster, cluster[1:])
    ]
    mean_gap = sum(gaps) / len(gaps) if gaps else 0.0
    start = min(peak.start for peak in cluster)
    end = max(peak.end for peak in cluster)
    score = len(cluster) * math.exp(-mean_gap / score_decay)
    missing_stats = [peak.label for peak in cluster if peak.label not in peak_stats]
    if missing_stats:
        raise ValueError(
            "Missing peak statistics for cluster peaks: " + ", ".join(missing_stats)
        )
    mean_log2foldchange = sum(
        peak_stats[peak.label]["log2FoldChange"] for peak in cluster
    ) / len(cluster)
    mean_fdr = sum(peak_stats[peak.label]["FDR"] for peak in cluster) / len(cluster)
    mean_of_mean_cpm = sum(peak_mean_cpm[peak.label] for peak in cluster) / len(cluster)
    mean_of_var_cpm = sum(peak_var_cpm[peak.label] for peak in cluster) / len(cluster)
    summary = {
        "chrom": cluster[0].chrom,
        "start": start,
        "end": end,
        "peaks": "/".join(peak.label for peak in cluster),
        "score": score,
        "peak_count": len(cluster),
        "cluster_span_bp": end - start,
        "mean_gap_bp": mean_gap,
        "mean_log2FoldChange": mean_log2foldchange,
        "mean_FDR": mean_fdr,
        "mean_of_mean_cpm": mean_of_mean_cpm,
        "mean_of_var_cpm": mean_of_var_cpm,
    }
    if any(peak.last_value is not None for peak in cluster):
        summary["Occurrence"] = "/".join(
            "" if peak.last_value is None else peak.last_value for peak in cluster
        )
    return summary


def format_cluster(summary):
    fields = [
        summary["chrom"],
        str(summary["start"]),
        str(summary["end"]),
        summary["peaks"],
        f'{summary["score"]:.6f}',
        str(summary["peak_count"]),
        str(summary["cluster_span_bp"]),
        f'{summary["mean_gap_bp"]:.2f}',
        f'{summary["mean_log2FoldChange"]:.6f}',
        f'{summary["mean_FDR"]:.6e}',
        f'{summary["mean_of_mean_cpm"]:.6f}',
        f'{summary["mean_of_var_cpm"]:.6f}',
    ]
    if "Occurrence" in summary:
        fields.append(summary["Occurrence"])
    return fields


def main():
    args = parse_args()
    if args.max_gap < 0:
        raise ValueError("--max-gap must be >= 0")
    if args.min_peaks < 1:
        raise ValueError("--min-peaks must be >= 1")
    score_decay = args.score_decay if args.score_decay is not None else args.max_gap
    if score_decay <= 0:
        raise ValueError("--score-decay must be > 0 (required when --max-gap is 0)")

    peaks = load_peaks(args.input)
    selected_peak_ids = {peak.label for peak in peaks}
    peak_stats = load_peak_stats(args.stats_tsv)
    peak_mean_cpm = load_peak_matrix_means(args.mean_cpm_tsv, selected_peak_ids)
    peak_var_cpm = load_peak_matrix_means(args.var_cpm_tsv, selected_peak_ids)
    clusters = [
        cluster
        for cluster in cluster_peaks(peaks, args.max_gap)
        if len(cluster) >= args.min_peaks
    ]
    summaries = [
        summarize_cluster(
            cluster, score_decay, peak_stats, peak_mean_cpm, peak_var_cpm
        )
        for cluster in clusters
    ]
    summaries.sort(
        key=lambda summary: (
            -summary["score"],
            -summary["peak_count"],
            chromosome_sort_key(summary["chrom"]),
            summary["start"],
            summary["end"],
        )
    )
    output_columns = [
        "chrom",
        "start",
        "end",
        "peaks",
        "score",
        "peak_count",
        "cluster_span_bp",
        "mean_gap_bp",
        "mean_log2FoldChange",
        "mean_FDR",
        "mean_of_mean_cpm",
        "mean_of_var_cpm",
    ]
    if any("Occurrence" in summary for summary in summaries):
        output_columns.append("Occurrence")

    with Path(args.output).open("w") as handle:
        handle.write(
            "\t".join(output_columns) + "\n"
        )
        for summary in summaries:
            handle.write("\t".join(format_cluster(summary)) + "\n")

    print(f"Read {len(peaks)} peaks from {args.input}")
    print(f"Wrote {len(clusters)} clusters to {args.output}")
    print(
        "Columns: "
        + ", ".join(output_columns[:-1] if output_columns[-1] == "Occurrence" else output_columns)
        + (", Occurrence" if output_columns[-1] == "Occurrence" else "")
    )


if __name__ == "__main__":
    main()
