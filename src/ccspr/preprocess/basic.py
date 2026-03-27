import numpy as np


def log_normalize_counts(x_counts: np.ndarray, scale: float = 1e4) -> np.ndarray:
    x = x_counts.astype(np.float64)
    lib = np.sum(x, axis=1, keepdims=True)
    lib[lib <= 0] = 1.0
    x = (x / lib) * float(scale)
    return np.log1p(x)


def select_top_variable_genes(
    x: np.ndarray,
    top_n: int = 5000,
) -> tuple[np.ndarray, np.ndarray]:
    n_features = x.shape[1]
    k = min(int(top_n), n_features)
    variances = np.var(x, axis=0)
    idx = np.argsort(-variances)[:k]
    return x[:, idx], idx


def standardize(
    x_train: np.ndarray,
    x_test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mu = np.mean(x_train, axis=0)
    sigma = np.std(x_train, axis=0)
    sigma[sigma < 1e-8] = 1.0
    return (x_train - mu) / sigma, (x_test - mu) / sigma
