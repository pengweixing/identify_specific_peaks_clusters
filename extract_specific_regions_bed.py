#!/usr/bin/env python3
"""Convert a subtype-specific differential TSV to a three-column BED file."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--specific-tsv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_peak_id(peak_id: str) -> tuple[str, int, int]:
    try:
        chrom, coordinates = peak_id.split(":", 1)
        start_text, end_text = coordinates.split("-", 1)
        start, end = int(start_text), int(end_text)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid peak_id: {peak_id}") from exc
    if not chrom or start < 0 or end <= start:
        raise ValueError(f"Invalid peak_id: {peak_id}")
    return chrom, start, end


def main() -> None:
    args = parse_args()
    seen: set[str] = set()
    count = 0
    with args.specific_tsv.open(newline="") as input_handle, args.output.open("w") as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        required = {"peak_id", "target_group"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{args.specific_tsv}: missing columns: {sorted(missing)}")
        for row in reader:
            peak_id = row["peak_id"]
            if peak_id in seen:
                raise ValueError(f"{args.specific_tsv}: duplicate peak_id: {peak_id}")
            seen.add(peak_id)
            chrom, start, end = parse_peak_id(peak_id)
            output_handle.write(f"{chrom}\t{start}\t{end}\n")
            count += 1
    print(f"Wrote {count:,} specific regions to {args.output}")


if __name__ == "__main__":
    main()
