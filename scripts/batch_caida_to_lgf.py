#!/usr/bin/env python3
"""
Batch convert all CAIDA AS relationship files to LGF format.
"""

import subprocess
import sys
from pathlib import Path
import os

def main():
    base_dir = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting")
    input_dir = base_dir / "experiments/datasets/caida/as-caida"
    output_dir = base_dir / "experiments/datasets/caida_lgf"
    script = base_dir / "scripts/caida_to_lgf.py"

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all .txt files
    txt_files = sorted(input_dir.glob("*.txt"))

    print(f"Found {len(txt_files)} CAIDA files to convert")
    print(f"Output directory: {output_dir}\n")

    success_count = 0
    fail_count = 0

    for i, input_file in enumerate(txt_files, 1):
        filename = input_file.stem
        output_file = output_dir / f"{filename}.lgf"

        print(f"[{i}/{len(txt_files)}] Converting {filename}...")

        try:
            result = subprocess.run(
                ["python3", str(script), "-i", str(input_file), "-o", str(output_file)],
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode == 0:
                print(f"  ✓ Success: {output_file.name}")
                success_count += 1
            else:
                print(f"  ✗ Failed with error:")
                print(f"    {result.stderr}")
                fail_count += 1

        except subprocess.TimeoutExpired:
            print(f"  ✗ Timeout (>5 minutes)")
            fail_count += 1
        except Exception as e:
            print(f"  ✗ Exception: {e}")
            fail_count += 1

        print()

    print("=" * 60)
    print(f"Conversion Summary:")
    print(f"  Successful: {success_count}/{len(txt_files)}")
    print(f"  Failed:     {fail_count}/{len(txt_files)}")
    print(f"  Output dir: {output_dir}")
    print("=" * 60)

    return 0 if fail_count == 0 else 1

if __name__ == '__main__':
    sys.exit(main())

