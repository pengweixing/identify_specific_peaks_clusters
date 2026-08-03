#!/usr/bin/env python3

import argparse
from concurrent.futures import ThreadPoolExecutor
import math
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

DEFAULT_CLUSTER_FILE = (
    "brca_subtype_vs_all_rest_deseq2/"
    "Basal_specific_regions.distance_clusters.ranked.bed"
)


def natural_key(path):
    name = Path(path).name
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", name)]


def collect_bigwigs(brca_dir, other_dir, max_other_tracks):
    brca_paths = sorted(Path(brca_dir).glob("*.bw"), key=natural_key)
    other_paths = sorted(Path(other_dir).glob("*.bw"), key=natural_key)
    if max_other_tracks > 0:
        other_paths = other_paths[:max_other_tracks]

    if not brca_paths:
        raise FileNotFoundError(f"No .bw files found in {brca_dir}")
    if not other_paths:
        raise FileNotFoundError(f"No .bw files found in {other_dir}")

    def group_name(path):
        return path.name.split("_", 1)[0]

    tracks = []
    for path in brca_paths + other_paths:
        tracks.append(
            {
                "path": path,
                "label": path.stem,
                "group": group_name(path),
                "is_brca": path in brca_paths,
            }
        )
    return tracks


def matching_chrom_from_chroms(chroms, chrom):
    if chrom in chroms:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in chroms:
        return chrom[3:]
    alt = f"chr{chrom}"
    if alt in chroms:
        return alt
    return None


def read_track_values(path, chrom, start, end, bins):
    with pyBigWig.open(str(path)) as bw:
        chroms = bw.chroms()
        bw_chrom = matching_chrom_from_chroms(chroms, chrom)
        if bw_chrom is None:
            return np.zeros(bins, dtype=float)
        chrom_len = chroms[bw_chrom]
        clipped_start = min(max(0, start), chrom_len)
        clipped_end = min(max(clipped_start + 1, end), chrom_len)
        values = bw.stats(
            bw_chrom,
            clipped_start,
            clipped_end,
            nBins=bins,
            type="mean",
        )
    return np.array(
        [0.0 if v is None or math.isnan(v) else float(v) for v in values],
        dtype=float,
    )


def parse_gtf_attrs(attrs):
    parsed = {}
    for item in attrs.rstrip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def load_gene_annotation(path):
    if not path:
        return {}
    annotation = Path(path)
    if not annotation.exists():
        raise FileNotFoundError(annotation)

    genes_by_chrom = {}
    gtf_gene_spans = {}
    with annotation.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue

            if len(fields) >= 9:
                if fields[2] not in {"gene", "transcript"}:
                    continue
                attrs = parse_gtf_attrs(fields[8])
                name = attrs.get("gene_name") or attrs.get("gene_id") or "."
                chrom = fields[0]
                start = int(fields[3]) - 1
                end = int(fields[4])
                strand = fields[6]
                gene_id = attrs.get("gene_id") or name
                gene_key = (chrom, gene_id, name, strand)
                if gene_key not in gtf_gene_spans:
                    gtf_gene_spans[gene_key] = {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "name": name,
                        "strand": strand,
                    }
                else:
                    gtf_gene_spans[gene_key]["start"] = min(
                        gtf_gene_spans[gene_key]["start"], start
                    )
                    gtf_gene_spans[gene_key]["end"] = max(
                        gtf_gene_spans[gene_key]["end"], end
                    )
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 and fields[3] else "."
            strand = fields[5] if len(fields) > 5 else "."
            genes_by_chrom.setdefault(chrom, []).append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "name": name,
                    "strand": strand,
                }
            )

    for gene in gtf_gene_spans.values():
        genes_by_chrom.setdefault(gene["chrom"], []).append(gene)
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(
            key=lambda item: (item["start"], item["end"], item["name"])
        )
    return genes_by_chrom


def overlapping_genes(genes_by_chrom, chrom, start, end):
    candidates = genes_by_chrom.get(chrom, [])
    if not candidates and chrom.startswith("chr"):
        candidates = genes_by_chrom.get(chrom[3:], [])
    if not candidates and not chrom.startswith("chr"):
        candidates = genes_by_chrom.get(f"chr{chrom}", [])
    return [gene for gene in candidates if gene["start"] < end and gene["end"] > start]


def assign_gene_lanes(genes):
    lane_ends = []
    assigned = []
    for gene in sorted(genes, key=lambda item: (item["start"], item["end"])):
        label_padding = max(1500, len(gene["name"]) * 180)
        padded_end = gene["end"] + label_padding
        lane = 0
        while lane < len(lane_ends) and gene["start"] <= lane_ends[lane]:
            lane += 1
        if lane == len(lane_ends):
            lane_ends.append(padded_end)
        else:
            lane_ends[lane] = padded_end
        assigned.append((gene, lane))
    return assigned, max(1, len(lane_ends))


def parse_peak_id(peak_id):
    chrom, coords = peak_id.split(":", 1)
    start, end = coords.split("-", 1)
    return chrom, int(start), int(end)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot IGV-like tracks for ranked Basal-specific peak clusters. "
            "Only peaks within each cluster are highlighted; CoRE regions are "
            "not plotted."
        )
    )
    parser.add_argument(
        "--cluster-file",
        default=DEFAULT_CLUSTER_FILE,
        help="Ranked cluster BED/TSV containing chrom, start, end, and peaks.",
    )
    parser.add_argument(
        "--brca-bigwig-dir",
        default="BRCA_subtype_random5_bigwigs",
        help="Directory containing BRCA subtype bigWig files.",
    )
    parser.add_argument(
        "--other-bigwig-dir",
        default="selected_10",
        help="Directory containing non-BRCA selected bigWig files.",
    )
    parser.add_argument(
        "--outdir",
        default="Basal_specific_cluster_igv_plots",
        help="Output directory for PNG figures.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top-ranked clusters to plot. Use 0 for all clusters.",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=10000,
        help="Flanking region drawn on each side of the cluster, in bp.",
    )
    parser.add_argument(
        "--scale-flank",
        type=int,
        default=0,
        help=(
            "Extra flank around the cluster used for y-axis scaling. "
            "Use 0 to scale on the cluster interval only."
        ),
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=1200,
        help="Number of bigWig summary bins across the displayed region.",
    )
    parser.add_argument(
        "--scale-bins",
        type=int,
        default=600,
        help="Number of bigWig summary bins used for scaling.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Parallel bigWig readers per figure.",
    )
    parser.add_argument(
        "--scale-mode",
        choices=("target", "all", "track"),
        default="target",
        help=(
            "Deprecated compatibility option. Plotting now uses hybrid scaling: "
            "each Basal track auto-scales to its own local maximum, while non-Basal "
            "tracks share a fixed scale equal to the average Basal maximum."
        ),
    )
    parser.add_argument(
        "--max-other-tracks",
        type=int,
        default=0,
        help="Limit non-BRCA tracks for quick previews. Use 0 for all tracks.",
    )
    parser.add_argument(
        "--gene-annotation",
        default="hg38.refGene.gtf",
        help="GTF or BED-like gene annotation. Use an empty string to skip.",
    )
    parser.add_argument(
        "--track-height",
        type=float,
        default=0.26,
        help="Figure height in inches per bigWig track.",
    )
    parser.add_argument(
        "--gene-track-height",
        type=float,
        default=1.7,
        help="Extra figure height in inches reserved for the gene track.",
    )
    parser.add_argument(
        "--signal-scale",
        type=float,
        default=1.5,
        help="Multiplier for signal height within each track lane.",
    )
    parser.add_argument("--fig-width", type=float, default=12.0)
    parser.add_argument(
        "--fig-height",
        type=float,
        default=0.0,
        help="Figure height in inches. Use 0 to calculate automatically.",
    )
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument(
        "--hide-title",
        action="store_true",
        help="Do not draw the top title text.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate and summarize the cluster file without reading bigWigs.",
    )
    return parser.parse_args()


def parse_cluster_peaks(value, expected_chrom):
    peaks = []
    for peak_id in str(value).split("/"):
        chrom, start, end = parse_peak_id(peak_id)
        if chrom != expected_chrom:
            raise ValueError(
                f"Cluster on {expected_chrom} contains peak on {chrom}: {peak_id}"
            )
        if start < 0 or end <= start:
            raise ValueError(f"Invalid peak interval: {peak_id}")
        peaks.append((chrom, start, end))
    if not peaks:
        raise ValueError("Cluster does not contain any peaks")
    return peaks


def load_clusters(path, top_n):
    clusters = pd.read_csv(path, sep="\t")
    required = {"chrom", "start", "end", "peaks"}
    missing = required - set(clusters.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")
    if "basal_specific_rank" in clusters.columns:
        clusters = clusters.sort_values("basal_specific_rank", kind="mergesort")
    if top_n > 0:
        clusters = clusters.head(top_n)

    parsed_peaks = []
    for row in clusters.itertuples(index=False):
        if int(row.start) < 0 or int(row.end) <= int(row.start):
            raise ValueError(f"Invalid cluster interval: {row.chrom}:{row.start}-{row.end}")
        peaks = parse_cluster_peaks(row.peaks, row.chrom)
        if any(start < row.start or end > row.end for _, start, end in peaks):
            raise ValueError(
                f"Cluster peaks extend outside interval: {row.chrom}:{row.start}-{row.end}"
            )
        parsed_peaks.append(peaks)
    clusters = clusters.copy()
    clusters["_parsed_peaks"] = parsed_peaks
    return clusters


def ordered_tracks_for_basal(tracks):
    target_tracks = [track for track in tracks if track["group"] == "Basal"]
    if not target_tracks:
        raise ValueError("No Basal bigWig tracks found")
    return target_tracks + [track for track in tracks if track["group"] != "Basal"]


def compute_track_denominators(ordered_tracks, scale_values_by_track):
    basal_maxima = [
        float(np.nanmax(values))
        for track, values in zip(ordered_tracks, scale_values_by_track)
        if track["group"] == "Basal" and values.size
    ]
    fixed_other_denominator = max(
        (sum(basal_maxima) / len(basal_maxima)) if basal_maxima else 0.0,
        1e-9,
    )
    denominators = []
    for track, values in zip(ordered_tracks, scale_values_by_track):
        if track["group"] == "Basal":
            ymax = float(np.nanmax(values)) if values.size else 0.0
            denominators.append(max(ymax, 1e-9))
        else:
            denominators.append(fixed_other_denominator)
    return denominators


def plot_genes(ax, genes_with_lanes, gene_lanes, start, end, window):
    gene_color = "#2f5d8c"
    for gene, lane in genes_with_lanes:
        y = gene_lanes - lane - 0.45
        clipped_start = max(start, gene["start"])
        clipped_end = min(end, gene["end"])
        if clipped_end <= clipped_start:
            continue
        ax.add_patch(
            plt.Rectangle(
                (clipped_start, y - 0.11),
                max(1, clipped_end - clipped_start),
                0.22,
                color=gene_color,
                alpha=0.9,
                linewidth=0,
            )
        )
        ax.text(
            clipped_start,
            y + 0.22,
            gene["name"],
            ha="left",
            va="bottom",
            fontsize=6,
            color="#1e3f63",
            clip_on=True,
        )
    ax.text(
        start - window * 0.006,
        gene_lanes * 0.5,
        "Genes",
        ha="right",
        va="center",
        fontsize=6,
        color="#1e3f63",
        fontweight="bold",
    )


def plot_cluster(row, tracks, genes_by_chrom, args, figure_rank):
    chrom = row["chrom"]
    cluster_start = int(row["start"])
    cluster_end = int(row["end"])
    peaks = row["_parsed_peaks"]
    start = max(0, cluster_start - args.flank)
    end = cluster_end + args.flank
    scale_start = max(0, cluster_start - args.scale_flank)
    scale_end = cluster_end + args.scale_flank
    window = end - start
    x = np.linspace(start, end, args.bins)

    ordered_tracks = ordered_tracks_for_basal(tracks)
    genes = overlapping_genes(genes_by_chrom, chrom, start, end) if genes_by_chrom else []
    genes_with_lanes, gene_lanes = assign_gene_lanes(genes) if genes else ([], 0)
    gene_space = gene_lanes + 1.0 if genes else 0.0
    signal_height = 0.82 * args.signal_scale
    lane_step = max(1.0, signal_height + 0.18)
    n_tracks = len(ordered_tracks)
    auto_fig_height = max(
        5.0,
        1.8
        + args.track_height * lane_step * n_tracks
        + (args.gene_track_height if genes else 0),
    )
    fig_height = args.fig_height if args.fig_height > 0 else auto_fig_height
    fig, ax = plt.subplots(figsize=(args.fig_width, fig_height))

    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as executor:
        values_by_track = list(
            executor.map(
                lambda track: read_track_values(
                    track["path"], chrom, start, end, args.bins
                ),
                ordered_tracks,
            )
        )
        scale_values_by_track = list(
            executor.map(
                lambda track: read_track_values(
                    track["path"], chrom, scale_start, scale_end, args.scale_bins
                ),
                ordered_tracks,
            )
        )

    denominators = compute_track_denominators(ordered_tracks, scale_values_by_track)

    track_base = gene_space
    for idx, (track, values, denominator) in enumerate(
        zip(ordered_tracks, values_by_track, denominators)
    ):
        offset = track_base + (n_tracks - idx - 1) * lane_step
        is_basal = track["group"] == "Basal"
        color = "#d62728" if is_basal else "#7f8790"
        alpha = 0.92 if is_basal else 0.62
        scaled = np.clip(values / denominator, 0, 1) * signal_height
        ax.fill_between(x, offset, offset + scaled, color=color, alpha=alpha, linewidth=0)
        ax.plot(x, offset + scaled, color=color, linewidth=0.25, alpha=0.95)
        ax.hlines(offset, start, end, color="#e5e5e5", linewidth=0.25)
        ax.text(
            start - window * 0.006,
            offset + signal_height * 0.46,
            track["label"],
            ha="right",
            va="center",
            fontsize=6,
            color="#222222" if is_basal else "#555555",
            fontweight="bold" if is_basal else "normal",
        )

    if genes:
        plot_genes(ax, genes_with_lanes, gene_lanes, start, end, window)

    for _, peak_start, peak_end in peaks:
        ax.axvspan(
            peak_start,
            peak_end,
            color="#f2c94c",
            alpha=0.38,
            linewidth=0,
            zorder=1,
        )

    ymax = track_base + n_tracks * lane_step + signal_height + 0.2
    ax.text(
        cluster_start,
        ymax - 0.16,
        f"{len(peaks)} cluster peaks",
        ha="left",
        va="top",
        fontsize=7,
        color="#9a7400",
        fontweight="bold",
    )
    ax.set_ylim(-0.2, ymax)
    ax.set_xlim(start, end)
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#9a9a9a")
    ax.set_xlabel(f"{chrom} position (bp)")
    ax.ticklabel_format(style="plain", axis="x")

    cluster_rank = int(row.get("basal_specific_rank", figure_rank))
    if not args.hide_title:
        title = (
            f"Basal cluster rank {cluster_rank}: {chrom}:{cluster_start:,}-{cluster_end:,}"
            f" | peaks={len(peaks)}"
        )
        if "basal_specific_rank_score" in row:
            title += f" | rank_score={float(row['basal_specific_rank_score']):.3g}"
        if "score" in row:
            title += f" | cluster_score={float(row['score']):.3g}"
        ax.set_title(title, fontsize=10, loc="left")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / (
        f"{cluster_rank:03d}_Basal_cluster_{chrom}_{cluster_start}_{cluster_end}.png"
    )
    fig.subplots_adjust(left=0.12, right=0.995, top=0.965, bottom=0.04)
    fig.savefig(outpath, dpi=args.dpi)
    plt.close(fig)
    return outpath


def main():
    args = parse_args()
    if args.top_n < 0:
        raise ValueError("--top-n must be >= 0")
    if args.flank < 0 or args.scale_flank < 0:
        raise ValueError("--flank and --scale-flank must be >= 0")
    if args.bins < 1 or args.scale_bins < 1:
        raise ValueError("--bins and --scale-bins must be >= 1")

    clusters = load_clusters(args.cluster_file, args.top_n)
    print(f"Read {len(clusters)} clusters from {args.cluster_file}")
    print(f"Cluster peaks: {sum(len(peaks) for peaks in clusters['_parsed_peaks'])}")
    if args.validate_only:
        print("Validation complete.")
        return

    tracks = collect_bigwigs(
        args.brca_bigwig_dir, args.other_bigwig_dir, args.max_other_tracks
    )
    genes_by_chrom = load_gene_annotation(args.gene_annotation)
    written = []
    for figure_rank, (_, row) in enumerate(clusters.iterrows(), start=1):
        outpath = plot_cluster(row, tracks, genes_by_chrom, args, figure_rank)
        written.append(outpath)
        print(outpath, flush=True)
    print(f"Done. Wrote {len(written)} figures to {args.outdir}")


if __name__ == "__main__":
    main()
