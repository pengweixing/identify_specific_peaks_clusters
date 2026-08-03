#!/usr/bin/env python3

import argparse
from collections import OrderedDict
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a peak-by-sample count matrix to CPM and summarize each peak "
            "by cancer type. Writes peak-by-cancer-type mean CPM and variance CPM "
            "matrices."
        )
    )
    parser.add_argument(
        "--matrix",
        default="data_m.Pancancer_peak.specific.mean_replicates.txt",
        help="Comma-separated peak-by-sample count matrix.",
    )
    parser.add_argument(
        "--sample-groups",
        default="sample_cancer_groups.tsv",
        help="Tab-separated file containing sample_name and cancer_group.",
    )
    parser.add_argument(
        "--mean-output",
        default="data_m.Pancancer_peak.specific.mean_replicates.mean_cpm_by_cancer_type.tsv",
        help="Output TSV for mean CPM values.",
    )
    parser.add_argument(
        "--variance-output",
        default="data_m.Pancancer_peak.specific.mean_replicates.var_cpm_by_cancer_type.tsv",
        help="Output TSV for CPM variance values.",
    )
    parser.add_argument(
        "--variance-ddof",
        type=int,
        default=0,
        help=(
            "Delta degrees of freedom for variance calculation (default: 0, "
            "population variance). Use 1 for sample variance."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Number of peaks processed at a time (default: 5000).",
    )
    return parser.parse_args()


def load_sample_groups(path):
    groups = pd.read_csv(path, sep="\t", dtype=str)
    required = {"sample_name", "cancer_group"}
    missing = required - set(groups.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")
    if groups[list(required)].isna().any().any():
        raise ValueError(f"{path}: sample_name and cancer_group must not be empty")
    duplicated = groups.loc[groups["sample_name"].duplicated(), "sample_name"].tolist()
    if duplicated:
        raise ValueError(f"{path}: duplicate sample names: {duplicated}")

    samples_by_group = OrderedDict()
    for row in groups.itertuples(index=False):
        samples_by_group.setdefault(row.cancer_group, []).append(row.sample_name)
    return groups["sample_name"].tolist(), samples_by_group


def load_matrix_header(path):
    columns = pd.read_csv(path, sep=",", nrows=0).columns.tolist()
    if not columns:
        raise ValueError(f"{path}: matrix has no columns")
    if columns[0] != "peak_id":
        raise ValueError(f"{path}: first column must be peak_id")
    sample_columns = columns[1:]
    if len(sample_columns) != len(set(sample_columns)):
        raise ValueError(f"{path}: matrix contains duplicate sample columns")
    return sample_columns


def validate_sample_names(matrix_samples, mapped_samples):
    matrix_set = set(matrix_samples)
    mapped_set = set(mapped_samples)
    missing_from_mapping = sorted(matrix_set - mapped_set)
    missing_from_matrix = sorted(mapped_set - matrix_set)
    if missing_from_mapping or missing_from_matrix:
        messages = []
        if missing_from_mapping:
            messages.append(f"matrix samples missing from mapping: {missing_from_mapping}")
        if missing_from_matrix:
            messages.append(f"mapped samples missing from matrix: {missing_from_matrix}")
        raise ValueError("; ".join(messages))


def numeric_values(chunk, sample_columns, matrix_path):
    values = chunk[sample_columns].apply(pd.to_numeric, errors="coerce")
    if values.isna().any().any():
        bad_columns = values.columns[values.isna().any()].tolist()
        raise ValueError(f"{matrix_path}: non-numeric values found in: {bad_columns}")
    if (values < 0).any().any():
        bad_columns = values.columns[(values < 0).any()].tolist()
        raise ValueError(f"{matrix_path}: negative values found in: {bad_columns}")
    return values


def calculate_library_sizes(matrix_path, sample_columns, chunk_size):
    library_sizes = pd.Series(0.0, index=sample_columns)
    peak_count = 0
    for chunk in pd.read_csv(matrix_path, sep=",", chunksize=chunk_size):
        values = numeric_values(chunk, sample_columns, matrix_path)
        library_sizes = library_sizes.add(values.sum(axis=0), fill_value=0)
        peak_count += len(chunk)
    bad_samples = library_sizes[library_sizes <= 0].index.tolist()
    if bad_samples:
        raise ValueError(f"{matrix_path}: non-positive library sizes: {bad_samples}")
    return library_sizes, peak_count


def summarize_chunk(
    chunk,
    sample_columns,
    samples_by_group,
    library_sizes,
    variance_ddof,
    matrix_path,
):
    values = numeric_values(chunk, sample_columns, matrix_path)
    cpm = values.divide(library_sizes, axis=1) * 1_000_000.0
    mean_summary = pd.DataFrame({"peak_id": chunk["peak_id"]})
    variance_summary = pd.DataFrame({"peak_id": chunk["peak_id"]})
    for cancer_group, group_samples in samples_by_group.items():
        group_cpm = cpm[group_samples]
        mean_summary[cancer_group] = group_cpm.mean(axis=1)
        variance_summary[cancer_group] = group_cpm.var(axis=1, ddof=variance_ddof)
    return mean_summary, variance_summary


def main():
    args = parse_args()
    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be >= 1")
    if args.variance_ddof < 0:
        raise ValueError("--variance-ddof must be >= 0")

    mapped_samples, samples_by_group = load_sample_groups(args.sample_groups)
    matrix_samples = load_matrix_header(args.matrix)
    validate_sample_names(matrix_samples, mapped_samples)
    if args.variance_ddof and any(
        len(samples) <= args.variance_ddof for samples in samples_by_group.values()
    ):
        small_groups = {
            group: len(samples)
            for group, samples in samples_by_group.items()
            if len(samples) <= args.variance_ddof
        }
        raise ValueError(
            "Not enough samples for requested --variance-ddof in groups: "
            f"{small_groups}"
        )

    library_sizes, peak_count = calculate_library_sizes(
        args.matrix, matrix_samples, args.chunk_size
    )
    mean_output = Path(args.mean_output)
    variance_output = Path(args.variance_output)
    for chunk_number, chunk in enumerate(
        pd.read_csv(args.matrix, sep=",", chunksize=args.chunk_size)
    ):
        mean_summary, variance_summary = summarize_chunk(
            chunk,
            matrix_samples,
            samples_by_group,
            library_sizes,
            args.variance_ddof,
            args.matrix,
        )
        write_header = chunk_number == 0
        mode = "w" if write_header else "a"
        mean_summary.to_csv(
            mean_output,
            sep="\t",
            index=False,
            mode=mode,
            header=write_header,
            float_format="%.10g",
        )
        variance_summary.to_csv(
            variance_output,
            sep="\t",
            index=False,
            mode=mode,
            header=write_header,
            float_format="%.10g",
        )

    print(f"Read {peak_count} peaks and {len(matrix_samples)} samples from {args.matrix}")
    print(f"Summarized {len(samples_by_group)} cancer types using CPM")
    print(f"Variance ddof: {args.variance_ddof}")
    print(f"Wrote mean CPM matrix: {mean_output}")
    print(f"Wrote variance CPM matrix: {variance_output}")


if __name__ == "__main__":
    main()
