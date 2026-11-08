"""Тесты на свойствах конвейера src/pipeline.py (направление B, свойства П1–П8)."""
import numpy as np
import pytest

from pipeline import (cpm_log, cv_accuracy, fit_predict_proba, make_folds,
                      prepare_features, select_top_variable)
from synthetic import make_synthetic

SEEDS = range(5)


@pytest.mark.parametrize("seed", SEEDS)
def test_p1_scaling_a_sample_does_not_change_cpm(seed):
    """П1: умножение всех каунтов образца на константу (глубина секвенирования) не меняет CPM."""
    counts, _ = make_synthetic(seed)
    mult = np.random.default_rng(100 + seed).integers(2, 10, size=counts.shape[1])
    assert np.allclose(cpm_log(counts), cpm_log(counts * mult), atol=1e-9)


@pytest.mark.parametrize("seed", SEEDS)
def test_p2_cpm_conserves_a_million_per_sample(seed):
    """П2: сохранение «массы»: после обратного преобразования CPM каждого образца даёт 1 000 000."""
    counts, _ = make_synthetic(seed)
    cpm = 2 ** cpm_log(counts) - 1
    assert np.allclose(cpm.sum(axis=0), 1e6, rtol=1e-9)


@pytest.mark.parametrize("seed", SEEDS)
def test_p3_physical_bounds_and_zero(seed):
    """П3: log2(CPM+1) лежит в [0; log2(1e6+1)], нулевой каунт даёт ровно 0."""
    counts, _ = make_synthetic(seed)
    v = cpm_log(counts)
    assert np.isfinite(v).all()
    assert v.min() >= 0
    assert v.max() <= np.log2(1e6 + 1) + 1e-9
    assert (v[counts == 0] == 0).all()


@pytest.mark.parametrize("seed", SEEDS)
def test_p4_order_of_samples_and_genes_is_irrelevant(seed):
    """П4: перестановка образцов и генов не меняет значения и набор отобранных генов."""
    counts, _ = make_synthetic(seed)
    rng = np.random.default_rng(200 + seed)
    ps, pg = rng.permutation(counts.shape[1]), rng.permutation(counts.shape[0])
    base = cpm_log(counts)
    assert np.allclose(cpm_log(counts[:, ps]), base[:, ps])
    assert np.allclose(cpm_log(counts[pg, :]), base[pg, :])
    sel = set(select_top_variable(base, 50))
    assert set(select_top_variable(base[:, ps], 50)) == sel
    assert set(pg[select_top_variable(base[pg, :], 50)]) == sel


@pytest.mark.parametrize("seed", SEEDS)
def test_p5_folds_are_a_stratified_partition(seed):
    """П5: каждый образец попадает в test ровно один раз, train и test не пересекаются,
    в каждом test-фолде есть все классы."""
    _, y = make_synthetic(seed)
    folds = make_folds(y, 5, seed)
    all_test = np.concatenate([te for _, te in folds])
    assert sorted(all_test) == list(range(len(y)))
    for tr, te in folds:
        assert not set(tr) & set(te)
        assert set(tr) | set(te) == set(range(len(y)))
        assert set(y[te]) == set(y)


@pytest.mark.parametrize("seed", SEEDS)
def test_p6_probabilities_are_a_distribution(seed):
    """П6: вероятности классов лежат в [0; 1], в каждой строке суммируются в 1,
    столбцы соответствуют классам обучающей выборки."""
    counts, y = make_synthetic(seed)
    X, _ = prepare_features(counts, k=100)
    tr, te = make_folds(y, 5, seed)[0]
    classes, proba = fit_predict_proba(X[tr], y[tr], X[te], seed=seed)
    assert proba.shape == (len(te), len(set(y[tr])))
    assert set(classes) == set(y[tr])
    assert (proba >= 0).all() and (proba <= 1).all()
    assert np.allclose(proba.sum(axis=1), 1.0)


def test_p7_same_seed_gives_same_result():
    """П7: детерминированность: тот же seed — то же предсказание и та же точность."""
    counts, y = make_synthetic(0)
    X, _ = prepare_features(counts, k=100)
    tr, te = make_folds(y, 5, 0)[0]
    _, p1 = fit_predict_proba(X[tr], y[tr], X[te], seed=3)
    _, p2 = fit_predict_proba(X[tr], y[tr], X[te], seed=3)
    assert np.array_equal(p1, p2)
    assert cv_accuracy(X, y, seed=3) == cv_accuracy(X, y, seed=3)


@pytest.mark.parametrize("seed", SEEDS)
def test_p8_more_reads_of_a_gene_never_lower_its_cpm(seed):
    """П8: монотонность: рост каунта гена в образце (при прочих равных) не уменьшает его CPM."""
    counts, _ = make_synthetic(seed)
    rng = np.random.default_rng(300 + seed)
    g, s = rng.integers(counts.shape[0]), rng.integers(counts.shape[1])
    more = counts.copy()
    more[g, s] += 5
    assert cpm_log(more)[g, s] >= cpm_log(counts)[g, s]


def test_p9_negative_control_labels_matter():
    """П9: отрицательный контроль: на настоящих метках точность высокая, на перемешанных —
    около уровня случая. Конвейер не «угадывает» метки по артефактам."""
    counts, y = make_synthetic(0, n_per_class=40)
    X, _ = prepare_features(counts, k=100)
    y_shuf = np.random.default_rng(0).permutation(y)
    assert cv_accuracy(X, y, seed=0) >= 0.95
    assert cv_accuracy(X, y_shuf, seed=0) <= 0.55
