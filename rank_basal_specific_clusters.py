#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_INPUT = "brca_subtype_vs_all_rest_deseq2/Basal2.distance_clusters.bed"
DEFAULT_OUTPUT = (
    "brca_subtype_vs_all_rest_deseq2/"
    "Basal2.distance_clusters.basal_specific_ranked.bed"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rank Basal peak clusters using a weighted percentile score. The "
            "structural cluster score receives the largest default weight."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input cluster BED.")
    parser.add_argument(
        "--output", default=DEFAULT_OUTPUT, help="Output ranked cluster BED."
    )
    parser.add_argument(
        "--score-weight",
        type=float,
        default=0.5,
        help="Weight for the structural cluster score percentile (default: 0.40).",
    )
    parser.add_argument(
        "--log2fc-weight",
        type=float,
        default=0.15,
        help="Weight for mean log2 fold-change percentile (default: 0.15).",
    )
    parser.add_argument(
        "--fdr-weight",
        type=float,
        default=0.15,
        help="Weight for -log10(mean FDR) percentile (default: 0.15).",
    )
    parser.add_argument(
        "--mean-cpm-weight",
        type=float,
        default=0.1,
        help="Weight for mean CPM percentile (default: 0.25).",
    )
    parser.add_argument(
        "--var-cpm-weight",
        type=float,
        default=0.05,
        help="Weight for CPM variance percentile (default: 0.05).",
    )
    parser.add_argument(
        "--variance-preference",
        choices=("low", "high"),
        default="low",
        help=(
            "Whether lower or higher mean CPM variance receives a better rank "
            "(default: low)."
        ),
    )
    return parser.parse_args()


def validate_input(df, path):
    required_columns = {
        "score",
        "mean_log2FoldChange",
        "mean_FDR",
        "mean_of_mean_cpm",
        "mean_of_var_cpm",
    }
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"{path}: missing required columns: {sorted(missing_columns)}")

    for column in required_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")
        if not np.isfinite(df[column]).all():
            raise ValueError(f"{path}: {column} contains missing or non-finite values")
    if (df["mean_FDR"] < 0).any():
        raise ValueError(f"{path}: mean_FDR must be >= 0")


def percentile_rank(values, higher_is_better=True):
    return values.rank(pct=True, method="average", ascending=higher_is_better)


def main():
    args = parse_args()
    weights = {
        "score_percentile": args.score_weight,
        "mean_log2FoldChange_percentile": args.log2fc_weight,
        "fdr_significance_percentile": args.fdr_weight,
        "mean_of_mean_cpm_percentile": args.mean_cpm_weight,
        "mean_of_var_cpm_stability_percentile": args.var_cpm_weight,
    }
    if any(weight < 0 for weight in weights.values()):
        raise ValueError("All weights must be >= 0")
    weight_sum = sum(weights.values())
    if weight_sum <= 0:
        raise ValueError("At least one weight must be > 0")
    weights = {name: weight / weight_sum for name, weight in weights.items()}

    df = pd.read_csv(args.input, sep="\t")
    if df.empty:
        raise ValueError(f"{args.input}: no clusters found")
    validate_input(df, args.input)

    df["neg_log10_mean_FDR"] = -np.log10(df["mean_FDR"].clip(lower=1e-300))
    df["score_percentile"] = percentile_rank(df["score"])
    df["mean_log2FoldChange_percentile"] = percentile_rank(
        df["mean_log2FoldChange"]
    )
    df["fdr_significance_percentile"] = percentile_rank(df["neg_log10_mean_FDR"])
    df["mean_of_mean_cpm_percentile"] = percentile_rank(df["mean_of_mean_cpm"])
    df["mean_of_var_cpm_stability_percentile"] = percentile_rank(
        df["mean_of_var_cpm"],
        higher_is_better=args.variance_preference == "high",
    )

    df["basal_specific_rank_score"] = sum(
        df[column] * weight for column, weight in weights.items()
    )
    df = df.sort_values(
        by=[
            "basal_specific_rank_score",
            "score",
            "mean_log2FoldChange",
            "neg_log10_mean_FDR",
            "chrom",
            "start",
            "end",
        ],
        ascending=[False, False, False, False, True, True, True],
        kind="mergesort",
    )
    df.insert(0, "basal_specific_rank", range(1, len(df) + 1))
    df.to_csv(args.output, sep="\t", index=False, float_format="%.10g")

    print(f"Read {len(df)} clusters from {args.input}")
    print(f"Wrote ranked clusters: {Path(args.output)}")
    print("Normalized weights:")
    for name, weight in weights.items():
        print(f"  {name}: {weight:.3f}")
    print(f"Variance preference: {args.variance_preference}")


if __name__ == "__main__":
    main()
