
"""Backward-compatible CLI for the E-Routing scenario runtime."""

from runtime import ScenarioRequest, TrustedRun, execute, main

__all__ = ["ScenarioRequest", "TrustedRun", "execute", "main"]


if __name__ == "__main__":
    main()
