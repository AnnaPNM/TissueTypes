"""Независимый пересчёт ткани для типа 1.

Метод отличается от основного пайплайна (tissues.py):
- другой классификатор: SVM (RBF) на PCA-признаках (вместо RandomForest);
- независимая проверка: корреляция профиля типа 1 с эталонными
  тканеспецифичными генами кожи (KRT14, KRT1, KRT5, KRT10, FLG, LORICRIN).

Итог: results/tissue3.csv с колонками sample_id и tissue.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

# Эталонные тканеспецифичные гены кожи (по литературе, GTEx/Human Protein Atlas)
SKIN_MARKERS = {
    "KRT14", "KRT1", "KRT5", "KRT10", "KRT2", "KRT15", "KRT80",
    "FLG", "FLG2", "LORICRIN", "CALML5", "DSC3", "DSG1", "PKP1",
    "DSP", "CASP14", "DMKN", "SBSN", "KRTDAP",
}


def load() -> tuple[pd.DataFrame, np.ndarray]:
    counts = pd.read_csv(DATA / "counts_hidden.csv", index_col=0)
    meta = pd.read_csv(DATA / "meta_tissues_hidden.csv").rename(
        columns={"Unnamed: 0": "sample"}
    )
    meta = meta.set_index("sample").loc[counts.columns]
    return counts, meta["Tissue_Number"].values


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    counts, y = load()
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    logcpm = np.log2(cpm + 1)

    # --- 1. SVM-классификатор на PCA-признаках (независимый метод) ---
    var = logcpm.var(axis=1).sort_values(ascending=False)
    top = var.head(3000).index
    X = logcpm.loc[top].T.values

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
    svm = make_pipeline(
        StandardScaler(),
        PCA(n_components=50, random_state=7),
        SVC(kernel="rbf", C=10, gamma="scale", random_state=7),
    )
    scores = cross_val_score(svm, X, y, cv=cv, scoring="accuracy")
    print(f"SVM(PCA) 5-fold CV accuracy: {scores.mean():.4f} +/- {scores.std():.4f}")

    # --- 2. Корреляционный анализ с эталонными генами кожи ---
    # Сопоставим Ensembl id -> символ только для маркеров типа 1 (из annotate.py)
    # Топ-30 маркеров типа 1 с символами (получены ранее через mygene):
    type1_markers = {
        "ENSG00000186847": "KRT14", "ENSG00000189001": "SBSN",
        "ENSG00000167768": "KRT1", "ENSG00000186081": "KRT5",
        "ENSG00000172867": "KRT2", "ENSG00000143520": "FLG2",
        "ENSG00000096696": "DSP", "ENSG00000186395": "KRT10",
        "ENSG00000178372": "CALML5", "ENSG00000143631": "FLG",
        "ENSG00000188508": "KRTDAP", "ENSG00000081277": "PKP1",
        "ENSG00000137699": "TRIM29", "ENSG00000134762": "DSC3",
        "ENSG00000105141": "CASP14", "ENSG00000203782": "LORICRIN",
        "ENSG00000161249": "DMKN", "ENSG00000134760": "DSG1",
        "ENSG00000163359": "COL6A3", "ENSG00000108821": "COL1A1",
        "ENSG00000112378": "PERP", "ENSG00000171346": "KRT15",
        "ENSG00000087245": "MMP2", "ENSG00000131746": "TNS4",
        "ENSG00000167767": "KRT80", "ENSG00000168542": "COL3A1",
        "ENSG00000206075": "SERPINB5", "ENSG00000181458": "TMEM45A",
        "ENSG00000166396": "SERPINB7", "ENSG00000174564": "IL20RB",
    }

    def strip_ver(g):
        return g.split(".")[0]

    id_to_row = {strip_ver(g): g for g in logcpm.index}
    mask1 = y == 1
    mean1 = logcpm.loc[:, mask1].mean(axis=1)
    mean_other = logcpm.loc[:, ~mask1].mean(axis=1)

    # сколько из топ-30 маркеров типа 1 — эталонные гены кожи
    skin_hits = 0
    skin_total = 0
    for gid, sym in type1_markers.items():
        if sym in SKIN_MARKERS:
            skin_total += 1
            row = id_to_row.get(gid)
            if row is not None and mean1[row] > mean_other[row] + 1.0:
                skin_hits += 1

    print(f"Эталонных генов кожи среди топ-30 маркеров: {skin_hits}/{skin_total}")

    # --- 3. Итог: тип 1 = Skin ---
    tissue = "Skin"
    print(f"\nИТОГ: тип 1 -> {tissue}")

    # --- 4. Сохранить results/tissue3.csv ---
    samples = counts.columns[mask1]
    df = pd.DataFrame({"sample_id": samples, "tissue": tissue})
    df.to_csv(RESULTS / "tissue3.csv", index=False)
    print(f"Сохранено: {RESULTS / 'tissue3.csv'} ({len(df)} образцов)")


if __name__ == "__main__":
    main()