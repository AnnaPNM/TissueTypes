"""Маленькие синтетические данные с известной структурой для тестов на свойствах."""
import numpy as np


def make_synthetic(seed=0, n_genes=300, n_per_class=30, n_classes=3, markers_per_class=30):
    """Каунты (гены x образцы) и метки 1..n_classes.

    У каждого класса свой набор из markers_per_class генов с 8-кратным повышением;
    размеры библиотек образцов различаются в 4 раза.
    """
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(1, n_classes + 1), n_per_class)
    base = rng.gamma(2.0, 20.0, size=n_genes)
    mu = np.tile(base[:, None], (1, labels.size))
    for c in range(n_classes):
        rows = slice(c * markers_per_class, (c + 1) * markers_per_class)
        mu[rows, labels == c + 1] *= 8
    depth = rng.uniform(0.5, 2.0, size=labels.size)
    counts = rng.poisson(mu * depth).astype(np.int64)
    return counts, labels
