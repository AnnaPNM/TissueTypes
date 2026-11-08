"""Ядро конвейера определения типа ткани по RNA-seq (то, что проверяет checks/).

Цепочка: каунты (гены x образцы) -> log2(CPM+1) -> top-k самых вариабельных генов ->
StandardScaler -> PCA -> логистическая регрессия. Отбор генов не использует метки.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

N_TOP_GENES = 3000
N_PCS = 50

ID_NAMES = {"sample_id", "sample", "samples", "sampleid", "id", "sample id"}
LABEL_NAMES = {"tissue", "tissue_type", "type", "label", "class", "group", "tissue_id"}


def cpm_log(counts, dtype=np.float64):
    """log2(CPM+1). counts: гены (строки) x образцы (столбцы).

    CPM считается по сумме каунтов каждого ОБРАЗЦА (по столбцу, axis=0).
    """
    x = np.asarray(counts, dtype=dtype)
    lib = x.sum(axis=0, keepdims=True)
    return np.log2(x / lib * 1e6 + 1)


def select_top_variable(logexpr, k):
    """Индексы k генов (строк) с наибольшей дисперсией по образцам; метки не используются."""
    k = min(k, logexpr.shape[0])
    var = logexpr.var(axis=1)
    idx = np.argsort(-var, kind="stable")[:k]
    return np.sort(idx)


def prepare_features(counts, k=N_TOP_GENES, dtype=np.float64):
    """counts (гены x образцы) -> X (образцы x k генов), индексы выбранных генов."""
    le = cpm_log(counts, dtype=dtype)
    idx = select_top_variable(le, k)
    return le[idx].T, idx


def make_folds(y, n_splits=5, seed=0):
    """Стратифицированное разбиение на фолды: список (train_idx, test_idx)."""
    y = np.asarray(y)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(skf.split(np.zeros(len(y)), y))


def fit_predict_proba(X_train, y_train, X_test, seed=0, n_pcs=N_PCS):
    """Обучает модель на train, возвращает (классы, вероятности на test)."""
    n_pcs = int(min(n_pcs, X_train.shape[0] - 1, X_train.shape[1]))
    model = make_pipeline(
        StandardScaler(),
        PCA(n_components=n_pcs, random_state=seed),
        LogisticRegression(max_iter=1000),
    )
    model.fit(X_train, y_train)
    return model.classes_, model.predict_proba(X_test)


def cv_accuracy(X, y, seed=0, n_splits=5):
    """Точность по объединённым предсказаниям стратифицированной k-fold CV."""
    y = np.asarray(y)
    correct = 0
    for tr, te in make_folds(y, n_splits, seed):
        classes, proba = fit_predict_proba(X[tr], y[tr], X[te], seed=seed)
        correct += int((classes[proba.argmax(axis=1)] == y[te]).sum())
    return correct / len(y)


def load_counts(path):
    """Таблица каунтов: гены в строках, образцы в столбцах, первый столбец — id генов."""
    return pd.read_csv(path, index_col=0)


def load_meta(path):
    """Метаданные -> pd.Series (индекс — id образца, значения — номер ткани).

    Столбец с id — первый (или с именем из ID_NAMES), столбец с меткой — с именем из
    LABEL_NAMES, иначе последний. Если схема другая, поправьте ID_NAMES/LABEL_NAMES.
    """
    m = pd.read_csv(path)
    if m.shape[1] < 2:
        raise ValueError(f"{Path(path).name}: нужен столбец с id образцов и столбец с меткой")
    lower = {c: str(c).strip().lower() for c in m.columns}
    id_col = next((c for c in m.columns if lower[c] in ID_NAMES), m.columns[0])
    label_col = next((c for c in m.columns if lower[c] in LABEL_NAMES and c != id_col), None)
    if label_col is None:
        label_col = [c for c in m.columns if c != id_col][-1]
    s = m.set_index(id_col)[label_col]
    s.index = s.index.map(str)
    return s
