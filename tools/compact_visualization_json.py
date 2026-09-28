#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path


class CompactionError(RuntimeError):
    pass


def compact_visualization_json(
        source: Path,
        destination: Path,
) -> int:
    temporary_path = destination.with_suffix(
        destination.suffix + ".tmp"
    )

    skipping = False
    buffered_line: str | None = None
    removed = 0

    try:
        with source.open(
                "r",
                encoding="utf-8",
        ) as input_file, temporary_path.open(
            "w",
            encoding="utf-8",
        ) as output_file:
            for line in input_file:
                stripped = line.lstrip()

                if (
                        not skipping and
                        stripped.startswith('"commodities": [')
                ):
                    if buffered_line is None:
                        raise CompactionError(
                            "Missing property before commodities array."
                        )

                    newline = (
                        "\n"
                        if buffered_line.endswith("\n")
                        else ""
                    )

                    previous = (
                        buffered_line[:-1]
                        if newline
                        else buffered_line
                    ).rstrip()

                    if not previous.endswith(","):
                        raise CompactionError(
                            "Property before commodities has no comma."
                        )

                    output_file.write(
                        previous[:-1] + newline
                    )

                    buffered_line = None
                    skipping = True
                    removed += 1
                    continue

                if skipping:
                    if stripped.startswith("]"):
                        skipping = False

                    continue

                if buffered_line is not None:
                    output_file.write(buffered_line)

                buffered_line = line

            if skipping:
                raise CompactionError(
                    "Input ended inside a commodities array."
                )

            if buffered_line is not None:
                output_file.write(buffered_line)

        with temporary_path.open(
                "r",
                encoding="utf-8",
        ) as compacted_file:
            json.load(compacted_file)

        temporary_path.replace(destination)

        return removed

    except Exception:
        temporary_path.unlink(
            missing_ok=True
        )

        raise


def destination_for(source: Path) -> Path:
    return source.with_name(
        source.stem +
        ".compact" +
        source.suffix
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Remove per-commodity failure arrays from "
            "visualization JSON files."
        )
    )

    parser.add_argument(
        "files",
        nargs="+",
        type=Path,
    )

    arguments = parser.parse_args()

    for source in arguments.files:
        source = source.expanduser().resolve()

        if not source.is_file():
            raise SystemExit(
                f"File not found: {source}"
            )

        destination = destination_for(
            source
        )

        removed = compact_visualization_json(
            source,
            destination,
        )

        source_mb = (
                source.stat().st_size /
                (1024 * 1024)
        )

        destination_mb = (
                destination.stat().st_size /
                (1024 * 1024)
        )

        print(
            f"{source.name}: removed {removed} commodity arrays; "
            f"{source_mb:.1f} MB -> {destination_mb:.1f} MB"
        )

        print(
            f"Written: {destination}"
        )


if __name__ == "__main__":
    main()