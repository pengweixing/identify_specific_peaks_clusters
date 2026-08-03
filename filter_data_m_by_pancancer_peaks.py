#!/usr/bin/env python3
"""Filter a TCGA ATAC count matrix using peaks from a BED file.

The BED file is expected to contain at least three tab-delimited columns:
chromosome, start, and end. The matrix is expected to be CSV-formatted with
peak IDs in the first column using the form chr:start-end.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_peak_ids(bed_path: Path) -> set[str]:
    peak_ids: set[str] = set()

    with bed_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{bed_path}:{line_number}: expected at least 3 columns, "
                    f"found {len(fields)}"
                )

            chrom, start, end = fields[:3]
            peak_ids.add(f"{chrom}:{start}-{end}")

    return peak_ids


def filter_matrix(matrix_path: Path, output_path: Path, peak_ids: set[str]) -> int:
    matched_rows = 0

    with matrix_path.open("r", encoding="utf-8", newline="") as in_handle:
        with output_path.open("w", encoding="utf-8", newline="") as out_handle:
            reader = csv.reader(in_handle)
            writer = csv.writer(out_handle, lineterminator="\n")

            header = next(reader, None)
            if header is None:
                raise ValueError(f"{matrix_path} is empty")

            writer.writerow(header)

            for row in reader:
                if not row:
                    continue
                if row[0] in peak_ids:
                    writer.writerow(row)
                    matched_rows += 1

    return matched_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter data_m.txt rows to peaks listed in Pancancer_peak.specific.bed."
    )
    parser.add_argument(
        "--bed",
        default="Pancancer_peak.specific.bed",
        type=Path,
        help="Input BED file with peaks. Default: %(default)s",
    )
    parser.add_argument(
        "--matrix",
        default="data_m.txt",
        type=Path,
        help="Input CSV matrix with peak IDs in column 1. Default: %(default)s",
    )
    parser.add_argument(
        "--output",
        default="data_m.Pancancer_peak.specific.txt",
        type=Path,
        help="Filtered output matrix. Default: %(default)s",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    peak_ids = load_peak_ids(args.bed)
    matched_rows = filter_matrix(args.matrix, args.output, peak_ids)

    print(f"Loaded {len(peak_ids):,} peaks from {args.bed}")
    print(f"Wrote {matched_rows:,} matching rows to {args.output}")


if __name__ == "__main__":
    main()
