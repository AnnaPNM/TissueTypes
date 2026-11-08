"""Прогон с повторами: настоящие метки против перемешанных (отрицательный контроль).

python scripts/run_experiment.py --step dry|smoke|full

dry   — фрагмент, один seed: тесты проходят и конвейер различает типы;
smoke — фрагмент, пять seed: правило принятия выполнено;
full  — полные данные, пять seed, те же правила; допускается только после smoke.
Каждый шаг печатает «ПРОЙДЕН» или завершается кодом 1. Строки прогона дописываются
в results/runs_log.csv.
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pipeline import cv_accuracy, load_counts, load_meta, prepare_features  # noqa: E402

SEEDS = [0, 1, 2, 3, 4]   # записаны в METHOD.md до первого прогона
DELTA = 0.20              # практический порог разности точностей
SHUF_MARGIN = 0.10        # перемешанные метки: не выше доли крупнейшего класса + 0.10
ACC_MIN = {"dry": 0.90, "smoke": 0.90, "full": 0.99}
LOG = ROOT / "results" / "runs_log.csv"


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "no-git"


def run_tests(full):
    env = dict(os.environ, TISSUES_FULL="1" if full else "0")
    p = subprocess.run([sys.executable, "-m", "pytest", "checks", "-q"], cwd=ROOT, env=env)
    return p.returncode == 0


def load(full):
    if full:
        counts = load_counts(ROOT / "data" / "counts_hidden.csv")
        meta = load_meta(ROOT / "data" / "meta_tissues_hidden.csv")
    else:
        counts = load_counts(ROOT / "checks" / "fixtures" / "counts_fixture.csv")
        meta = load_meta(ROOT / "checks" / "fixtures" / "meta_fixture.csv")
    y = meta.loc[counts.columns.astype(str)].astype(int).to_numpy()
    return counts.to_numpy(), y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", choices=["dry", "smoke", "full"], required=True)
    step = ap.parse_args().step
    full = step == "full"

    print(f"== шаг {step}: тесты checks/ ({'полные данные' if full else 'фрагмент'}) ==")
    if not run_tests(full):
        print("шаг НЕ ПРОЙДЕН: тесты checks/ падают")
        sys.exit(1)

    counts, y = load(full)
    X, _ = prepare_features(counts, dtype=np.float32 if full else np.float64)
    majority = np.bincount(y).max() / len(y)
    seeds = SEEDS[:1] if step == "dry" else SEEDS
    LOG.parent.mkdir(exist_ok=True)
    new = not LOG.exists()
    rows = []
    with open(LOG, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["step", "seed", "variant", "accuracy", "n_samples", "seconds", "commit"])
        for s in seeds:
            y_shuf = np.random.default_rng(1000 + s).permutation(y)
            for variant, labels in (("true", y), ("shuffled", y_shuf)):
                t0 = time.time()
                acc = cv_accuracy(X, labels, seed=s)
                rows.append((s, variant, acc))
                w.writerow([step, s, variant, f"{acc:.6f}", len(y), f"{time.time() - t0:.1f}", git_commit()])
                print(f"seed={s} {variant:8s} accuracy={acc:.4f}")

    a_true = np.array([r[2] for r in rows if r[1] == "true"])
    a_shuf = np.array([r[2] for r in rows if r[1] == "shuffled"])
    d = a_true - a_shuf
    n = len(d)
    ok1 = a_true.mean() >= ACC_MIN[step]
    ok2 = a_shuf.mean() <= majority + SHUF_MARGIN
    if n >= 2:
        sd = d.std(ddof=1)
        half = stats.t.ppf(0.975, n - 1) * sd / np.sqrt(n)
        lo, hi = d.mean() - half, d.mean() + half
        ok3 = lo > DELTA
        print(f"\nнастоящие метки: среднее {a_true.mean():.4f}, sd {a_true.std(ddof=1):.4f}")
        print(f"перемешанные:    среднее {a_shuf.mean():.4f}, sd {a_shuf.std(ddof=1):.4f} "
              f"(крупнейший класс {majority:.4f})")
        print(f"парная разность: среднее {d.mean():.4f}, 95 %-интервал [{lo:.4f}; {hi:.4f}], порог {DELTA}")
    else:
        ok3 = d.mean() > DELTA
        print(f"\nнастоящие {a_true.mean():.4f}, перемешанные {a_shuf.mean():.4f} (один seed, интервала нет)")
    print(f"условие 1 (точность >= {ACC_MIN[step]}): {'да' if ok1 else 'НЕТ'}")
    print(f"условие 2 (перемешанные <= крупнейший класс + {SHUF_MARGIN}): {'да' if ok2 else 'НЕТ'}")
    print(f"условие 3 (нижняя граница разности > {DELTA}): {'да' if ok3 else 'НЕТ'}")
    if ok1 and ok2 and ok3:
        print(f"\nшаг {step}: ПРОЙДЕН")
        return
    print(f"\nшаг {step}: НЕ ПРОЙДЕН")
    sys.exit(1)


if __name__ == "__main__":
    main()
