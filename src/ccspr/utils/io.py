import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_yaml(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def save_json(obj: dict[str, Any], path: str | Path) -> None:
    p = Path(path)
    ensure_dir(p.parent)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def hash_array(x: np.ndarray, precision: int = 8) -> str:
    x = np.asarray(x)
    if x.dtype.kind == "f":
        x = np.round(x.astype(np.float64), precision)
    c = np.ascontiguousarray(x)
    h = hashlib.sha1()
    h.update(str(c.shape).encode("utf-8"))
    h.update(str(c.dtype).encode("utf-8"))
    h.update(c.tobytes())
    return h.hexdigest()
