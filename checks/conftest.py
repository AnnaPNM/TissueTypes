import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pipeline import load_counts, load_meta  # noqa: E402

FULL = os.environ.get("TISSUES_FULL") == "1"
FIXTURES = ROOT / "checks" / "fixtures"
EXPECTED_FULL_COUNTS = {1: 1809, 2: 578, 3: 2642, 4: 929, 5: 89}


@pytest.fixture(scope="session")
def data():
    """(каунты, метки, режим). По умолчанию — маленький фрагмент checks/fixtures/;
    TISSUES_FULL=1 — полные данные из data/."""
    if FULL:
        counts_p = ROOT / "data" / "counts_hidden.csv"
        meta_p = ROOT / "data" / "meta_tissues_hidden.csv"
        if not counts_p.exists() or not meta_p.exists():
            pytest.skip("полные данные data/ не найдены")
        mode = "full"
    else:
        counts_p = FIXTURES / "counts_fixture.csv"
        meta_p = FIXTURES / "meta_fixture.csv"
        if not counts_p.exists() or not meta_p.exists():
            pytest.fail("нет checks/fixtures/: запусти python checks/make_fixture.py")
        mode = "fixture"
    return load_counts(counts_p), load_meta(meta_p), mode


@pytest.fixture(scope="session")
def result_table():
    p = ROOT / "results" / "tissues_result.csv"
    if not p.exists():
        pytest.skip("results/tissues_result.csv не найден")
    import pandas as pd
    return pd.read_csv(p)
