#!/usr/bin/env bash
# Вносит ошибку в отдельной ветке и сохраняет вывод проверки «до» и «после».
# Запуск из корня проекта, все изменения закоммичены: bash scripts/catch_error.sh
set -u
cd "$(git rev-parse --show-toplevel)"
BASE=$(git rev-parse --abbrev-ref HEAD)
BUG=inject-orientation-bug
if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "Сначала закоммить изменения (git status)"; exit 2
fi
mkdir -p results
git checkout -q -b "$BUG" || exit 2
python3 - <<'PY'
import pathlib
p = pathlib.Path("src/pipeline.py")
s = p.read_text()
old = "lib = x.sum(axis=0, keepdims=True)"
new = "lib = x.sum(axis=1, keepdims=True)"
assert old in s, "в src/pipeline.py нет строки для подмены"
p.write_text(s.replace(old, new))
PY
git commit -qam "Внесённая ошибка: CPM по сумме гена (строки) вместо суммы образца (столбца)"
python3 -m pytest checks -q > /tmp/caught_error_before.txt 2>&1
echo "exit code: $?" >> /tmp/caught_error_before.txt
git checkout -q "$BASE"
python3 -m pytest checks -q > results/caught_error_after.txt 2>&1
echo "exit code: $?" >> results/caught_error_after.txt
cp /tmp/caught_error_before.txt results/caught_error_before.txt
echo "--- до (ветка $BUG) ---"; tail -n 12 results/caught_error_before.txt
echo "--- после (ветка $BASE) ---"; tail -n 4 results/caught_error_after.txt
