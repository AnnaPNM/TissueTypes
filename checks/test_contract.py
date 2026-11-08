"""Договор о данных: схема, типы, диапазоны, ключи, тождества, согласованность с результатом."""
import numpy as np
import pandas as pd

from conftest import EXPECTED_FULL_COUNTS

TYPES = {1, 2, 3, 4, 5}


def test_c1_counts_are_finite_nonnegative_integers(data):
    counts, _, _ = data
    v = counts.to_numpy()
    assert np.isfinite(v).all(), "в каунтах есть NaN или inf"
    assert (v >= 0).all(), "в каунтах есть отрицательные значения"
    assert np.all(v == np.floor(v)), "каунты должны быть целыми"


def test_c2_gene_and_sample_ids_are_unique_and_nonempty(data):
    counts, meta, _ = data
    for name, ids in (("гены", counts.index), ("образцы в каунтах", counts.columns),
                      ("образцы в метаданных", meta.index)):
        s = pd.Series(ids.astype(str))
        assert not s.duplicated().any(), f"повторяются id: {name}"
        assert (s.str.strip() != "").all(), f"пустые id: {name}"


def test_c3_every_sample_has_reads(data):
    counts, _, _ = data
    assert (counts.sum(axis=0) > 0).all(), "есть образец с нулевой суммой каунтов (CPM не определён)"


def test_c4_labels_are_integers_1_to_5(data):
    _, meta, _ = data
    assert meta.notna().all(), "есть пропущенные метки"
    v = meta.to_numpy()
    assert np.all(v == np.floor(v)), "метки должны быть целыми"
    assert set(int(x) for x in v) <= TYPES, "метки вне 1..5"


def test_c5_samples_match_between_counts_and_meta(data):
    counts, meta, _ = data
    assert set(counts.columns.astype(str)) == set(meta.index), "образцы в каунтах и метаданных не совпадают"


def test_c6_all_five_types_present_with_known_sizes(data):
    _, meta, mode = data
    sizes = meta.astype(int).value_counts().to_dict()
    assert set(sizes) == TYPES, "представлены не все пять типов"
    if mode == "full":
        assert sizes == EXPECTED_FULL_COUNTS, f"размеры типов {sizes} != {EXPECTED_FULL_COUNTS}"


def test_c7_result_table_agrees_with_meta(data, result_table):
    """Тождество: каждый зашифрованный тип сопоставлен ровно одной ткани, разным типам — разные ткани."""
    _, meta, mode = data
    r = result_table
    assert {"sample_id", "tissue"} <= set(r.columns)
    assert r["tissue"].notna().all() and (r["tissue"].astype(str).str.strip() != "").all()
    r = r.assign(sample_id=r["sample_id"].astype(str))
    assert not r["sample_id"].duplicated().any(), "повторяются sample_id в результате"
    missing = set(meta.index) - set(r["sample_id"])
    assert not missing, f"в результате нет {len(missing)} образцов из метаданных"
    if mode == "full":
        assert set(r["sample_id"]) == set(meta.index)
    j = r.set_index("sample_id").loc[meta.index].assign(type=meta.astype(int).to_numpy())
    per_type = j.groupby("type")["tissue"].nunique()
    assert (per_type == 1).all(), f"типу соответствует несколько тканей: {per_type.to_dict()}"
    assert j.groupby("type")["tissue"].first().nunique() == 5, "разным типам назначена одна ткань"
