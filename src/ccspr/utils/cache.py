import pickle
from pathlib import Path
from typing import Any

import numpy as np

from .io import ensure_dir


def cache_path(cache_dir: str | Path, prefix: str, key: str, suffix: str) -> Path:
    return ensure_dir(Path(cache_dir) / prefix) / f"{key}.{suffix}"


def load_npy(path: str | Path) -> np.ndarray | None:
    p = Path(path)
    if p.exists():
        return np.load(p)
    return None


def save_npy(path: str | Path, arr: np.ndarray) -> None:
    p = Path(path)
    ensure_dir(p.parent)
    np.save(p, arr)


def load_pickle(path: str | Path) -> Any:
    p = Path(path)
    if not p.exists():
        return None
    with open(p, "rb") as f:
        return pickle.load(f)


def save_pickle(path: str | Path, obj: Any) -> None:
    p = Path(path)
    ensure_dir(p.parent)
    with open(p, "wb") as f:
        pickle.dump(obj, f)
