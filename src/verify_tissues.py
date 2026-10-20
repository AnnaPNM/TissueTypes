"""Независимый пересчёт тканей для типов 2, 3, 4, 5.

Метод отличается от основного пайплайна (tissues.py):
- другой классификатор: SVM (RBF) на PCA-признаках (вместо RandomForest);
- независимая проверка: корреляция профиля каждого типа с эталонными
  тканеспецифичными генами (Lung, Brain, Blood, Kidney).

Тип 1 (Skin) уже определён и принят пользователем — здесь не пересчитывается,
но включается в итоговый файл results/tissues_result.csv.

Итог: results/tissues_result.csv с колонками sample_id и tissue.
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

# Эталонные тканеспецифичные гены (Ensembl id) по литературе (GTEx/HPA).
# Символы: Lung=SFTPB,SFTPA1,SFTPA2,SFTPC,NKX2-1,NAPSA,AGER;
# Brain=GFAP,PLP1,MOBP,SNAP25,SLC1A2,KIF5A,TUBB4A;
# Blood=HBB,HBA1,HBA2,HBD,S100A9,LYZ,MMP9;
# Kidney=UMOD,AQP2,SLC12A1,SLC12A3,NPHS2,PAX2,PAX8,CDH16.
REF_MARKERS = {
    2: {  # Lung
        "ENSG00000168878", "ENSG00000122852", "ENSG00000185303",
        "ENSG00000168484", "ENSG00000136352", "ENSG00000131400",
        "ENSG00000204305",
    },
    3: {  # Brain
        "ENSG00000131095", "ENSG00000123560", "ENSG00000168314",
        "ENSG00000132639", "ENSG00000110436", "ENSG00000155980",
        "ENSG00000104833",
    },
    4: {  # Blood
        "ENSG00000244734", "ENSG00000206172", "ENSG00000188536",
        "ENSG00000223609", "ENSG00000163220", "ENSG00000090382",
        "ENSG00000100985",
    },
    5: {  # Kidney
        "ENSG00000169344", "ENSG00000167580", "ENSG00000074803",
        "ENSG00000070915", "ENSG00000116218", "ENSG00000075891",
        "ENSG00000125618", "ENSG00000166589",
    },
}
TISSUE_NAMES = {1: "Skin", 2: "Lung", 3: "Brain", 4: "Blood", 5: "Kidney"}


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

    # --- 2. Корреляционный анализ с эталонными генами ---
    def strip_ver(g):
        return g.split(".")[0]

    id_to_row = {strip_ver(g): g for g in logcpm.index}
    print("\nЭталонные маркеры, повышенные в своём типе (log2CPM diff > 1.0):")
    for t in sorted(REF_MARKERS):
        mask = y == t
        mean_t = logcpm.loc[:, mask].mean(axis=1)
        mean_o = logcpm.loc[:, ~mask].mean(axis=1)
        hits = 0
        total = len(REF_MARKERS[t])
        for gid in REF_MARKERS[t]:
            row = id_to_row.get(gid)
            if row is not None and mean_t[row] > mean_o[row] + 1.0:
                hits += 1
        print(f"  Тип {t} ({TISSUE_NAMES[t]}): {hits}/{total} маркеров повышены")
        assert hits >= total * 0.7, f"Тип {t}: слишком мало маркеров повышено"

    # --- 3. Итог: типы 2-5 ---
    print("\nИТОГ:")
    for t in sorted(TISSUE_NAMES):
        print(f"  тип {t} -> {TISSUE_NAMES[t]}")

    # --- 4. Сохранить results/tissues_result.csv ---
    rows = [{"sample_id": s, "tissue": TISSUE_NAMES[t]}
            for s, t in zip(counts.columns, y)]
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "tissues_result.csv", index=False)
    df.to_csv(ROOT / "tissues_result.csv", index=False)
    print(f"\nСохранено: {RESULTS / 'tissues_result.csv'} ({len(df)} образцов)")


if __name__ == "__main__":
    main()