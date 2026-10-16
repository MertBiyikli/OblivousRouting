"""Compatibility CLI for canonical E-Routing analysis queries."""

from e_routing.analysis.queries import AnalysisAPI, main

__all__ = ["AnalysisAPI", "main"]


if __name__ == "__main__":
    main()