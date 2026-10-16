"""Compatibility CLI for E-Routing analysis validation."""

from e_routing.analysis.validation import main, near, validate_semantics

__all__ = ["main", "near", "validate_semantics"]


if __name__ == "__main__":
    main()