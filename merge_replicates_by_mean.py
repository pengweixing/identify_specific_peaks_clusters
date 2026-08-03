#!/usr/bin/env python3
"""Merge replicate columns in a TCGA ATAC matrix by mean value.

Input columns look like:
ACCx_025FE5F8_885E_433D_9018_7AE322A92285_X034_S09_L133_B1_T1_PMRG
ACCx_025FE5F8_885E_433D_9018_7AE322A92285_X034_S09_L134_B1_T2_PMRG

The merged sample barcode is the part before the run/lane fields:
ACCx_025FE5F8_885E_433D_9018_7AE322A92285
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def sample_barcode(column_name: str) -> str:
    fields = column_name.split("_")
    if len(fields) < 6:
        raise ValueError(f"Cannot parse sample barcode from column name: {column_name}")

    return "_".join(fields[:6])


def build_groups(header: list[str]) -> tuple[list[str], list[list[int]]]:
    barcodes: list[str] = []
    group_indices: list[list[int]] = []
    barcode_to_group: dict[str, int] = {}

    for index, column_name in enumerate(header):
        barcode = sample_barcode(column_name)
        group_index = barcode_to_group.get(barcode)

        if group_index is None:
            barcode_to_group[barcode] = len(barcodes)
            barcodes.append(barcode)
            group_indices.append([index])
        else:
            group_indices[group_index].append(index)

    return barcodes, group_indices


def mean_text(values: list[str]) -> str:
    numbers = [float(value) for value in values if value != ""]
    if not numbers:
        return ""

    mean_value = sum(numbers) / len(numbers)
    if mean_value.is_integer():
        return str(int(mean_value))

    return f"{mean_value:.6f}".rstrip("0").rstrip(".")


def merge_matrix(input_path: Path, output_path: Path) -> tuple[int, int, int]:
    row_count = 0

    with input_path.open("r", encoding="utf-8", newline="") as in_handle:
        with output_path.open("w", encoding="utf-8", newline="") as out_handle:
            reader = csv.reader(in_handle)
            writer = csv.writer(out_handle, lineterminator="\n")

            header = next(reader, None)
            if header is None:
                raise ValueError(f"{input_path} is empty")

            barcodes, group_indices = build_groups(header)
            writer.writerow(["peak_id", *barcodes])

            for row in reader:
                if not row:
                    continue

                peak_id = row[0]
                values = row[1:]
                merged_values = [
                    mean_text([values[index] for index in indices])
                    for indices in group_indices
                ]
                writer.writerow([peak_id, *merged_values])
                row_count += 1

    return len(header), len(barcodes), row_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge replicate matrix columns by sample barcode using mean values."
    )
    parser.add_argument(
        "--input",
        default="data_m.Pancancer_peak.specific.txt",
        type=Path,
        help="Input CSV matrix. Default: %(default)s",
    )
    parser.add_argument(
        "--output",
        default="data_m.Pancancer_peak.specific.mean_replicates.txt",
        type=Path,
        help="Output CSV matrix with merged replicate columns. Default: %(default)s",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_columns, output_columns, row_count = merge_matrix(args.input, args.output)

    print(f"Input sample columns: {input_columns:,}")
    print(f"Output merged sample columns: {output_columns:,}")
    print(f"Rows written: {row_count:,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
