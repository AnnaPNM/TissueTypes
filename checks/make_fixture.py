"""Делает маленький фрагмент данных для запуска верификатора без полных данных.

python checks/make_fixture.py [--per-type 20] [--genes 500] [--seed 0]
Читает data/counts_hidden.csv и data/meta_tissues_hidden.csv, пишет checks/fixtures/.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pipeline import load_meta  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-type", type=int, default=20)
    ap.add_argument("--genes", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    counts_p = ROOT / "data" / "counts_hidden.csv"
    meta = load_meta(ROOT / "data" / "meta_tissues_hidden.csv")
    header = pd.read_csv(counts_p, nrows=0)
    gene_col, samples = header.columns[0], [str(c) for c in header.columns[1:]]
    print(f"meta: {len(meta)} образцов, метки {sorted(meta.unique())}; в каунтах {len(samples)} образцов")
    print(f"столбец генов: {gene_col!r}; пример id образца: {samples[0]!r}; в метаданных: {meta.index[0]!r}")

    chosen = []
    for t in sorted(meta.unique()):
        ids = [s for s in meta.index[meta == t] if s in set(samples)]
        chosen += list(rng.choice(ids, size=min(a.per_type, len(ids)), replace=False))
    df = pd.read_csv(counts_p, index_col=0, usecols=[gene_col] + chosen)[chosen]

    logv = np.log2(df.to_numpy(dtype=float) + 1).var(axis=1)
    top = np.argsort(-logv)[: int(a.genes * 0.6)]
    rest = np.setdiff1d(np.arange(len(df)), top)
    rand = rng.choice(rest, size=a.genes - len(top), replace=False)
    df = df.iloc[np.sort(np.concatenate([top, rand]))]

    out = ROOT / "checks" / "fixtures"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "counts_fixture.csv")
    pd.DataFrame({"sample_id": chosen, "type": meta.loc[chosen].astype(int).to_numpy()}).to_csv(
        out / "meta_fixture.csv", index=False)
    print(f"фрагмент: {df.shape[0]} генов x {df.shape[1]} образцов -> {out}")


if __name__ == "__main__":
    main()
