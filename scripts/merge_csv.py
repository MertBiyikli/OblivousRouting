import argparse
from pathlib import Path

import pandas as pd


def merge_csv_files(input_dir: str | Path, output_file: str | Path) -> None:
    input_dir = Path(input_dir).expanduser().resolve()
    output_file = Path(output_file).expanduser().resolve()

    csv_files = sorted(input_dir.rglob("*.csv"))
    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return

    merged = pd.concat(
        (pd.read_csv(f) for f in csv_files),
        ignore_index=True,
    )
    merged.to_csv(output_file, index=False)
    print(f"Merged {len(csv_files)} CSV file(s) → {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge all CSV files in a directory into a single CSV."
    )
    parser.add_argument("input_dir", help="Directory to search for CSV files")
    parser.add_argument("output_file", help="Path for the merged output CSV")
    args = parser.parse_args()

    merge_csv_files(args.input_dir, args.output_file)
