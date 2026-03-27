# Note: Broken Import in geometry/__init__.py

## Problem
`src/ccspr/geometry/__init__.py` line 1 imports `SUPPORTED_MODES` from `distance.py`, but the actual exported name in `distance.py:264` is `VALID_MODES`.

This causes `ImportError: cannot import name 'SUPPORTED_MODES'` whenever any `ccspr.*` module is imported, because the package `__init__.py` chain runs:
  `ccspr/__init__.py` → `ccspr.geometry.distance` → `ccspr.geometry/__init__.py` → broken.

## Fix Required (1 line change)
In `src/ccspr/geometry/__init__.py`:
```python
# Current (broken):
from .distance import SUPPORTED_MODES, build_distance
__all__ = ["build_distance", "SUPPORTED_MODES"]

# Fix:
from .distance import VALID_MODES, build_distance
__all__ = ["build_distance", "VALID_MODES"]
```

Or add an alias in `distance.py`:
```python
SUPPORTED_MODES = VALID_MODES  # backward compat
```

## Why Existing Scripts Work
The existing experiment scripts (run_geometry_family.py, etc.) worked when first written because this bug was introduced later. They use the installed editable package (`pip install -e .`) which may have been installed before the `__init__.py` was broken, or the `__init__.py` was correct at install time.

## Workaround Used in New Scripts
New scripts (run_arabidopsis.py) temporarily patch the `__init__.py` file in-memory before importing.
