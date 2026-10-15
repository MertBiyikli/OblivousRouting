#!/usr/bin/env bash
set -euo pipefail

echo "Checking reproducible-build invariants..."

if [[ -f CMakeUserPresets.json ]]; then
  if git ls-files --error-unmatch CMakeUserPresets.json >/dev/null 2>&1; then
    echo "ERROR: CMakeUserPresets.json is tracked; it must remain machine-local." >&2
    exit 1
  fi
fi

# Every declared submodule must be initialized at the exact gitlink recorded by
# the parent repository.
if [[ -f .gitmodules ]]; then
  if git submodule status --recursive | grep -E '^[+-]' >/dev/null; then
    echo "ERROR: a submodule is missing or checked out at the wrong revision." >&2
    git submodule status --recursive >&2
    exit 1
  fi
fi

if [[ ! -f extern/amgcl/amgcl/backend/builtin.hpp ]]; then
  echo "ERROR: extern/amgcl is missing; initialize submodules recursively." >&2
  exit 1
fi

if [[ ! -f extern/xcut/CMakeLists.txt ]]; then
  echo "ERROR: extern/xcut is not present in a clean checkout." >&2
  echo "Track it as a pinned submodule or vendor it into the repository." >&2
  exit 1
fi

# Keep the project preset file machine-independent.
if grep -E '/opt/homebrew|/usr/local/lib/cmake|\$env\{HOME\}/local' CMakePresets.json >/dev/null; then
  echo "ERROR: machine-specific paths leaked into CMakePresets.json." >&2
  exit 1
fi

echo "Reproducibility checks passed."
