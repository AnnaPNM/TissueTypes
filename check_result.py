import csv, sys

path = "tissues_result.csv"
try:
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
except FileNotFoundError:
    sys.exit(f"FAIL: нет файла {path}")

if not rows:
    sys.exit("FAIL: файл пустой")
if not {"sample_id", "tissue"} <= set(rows[0].keys()):
    sys.exit("FAIL: нужны колонки sample_id и tissue")

tissues = [r["tissue"].strip() for r in rows]
if any(t == "" for t in tissues):
    sys.exit("FAIL: есть пустые значения tissue")
if len(set(tissues)) != 5:
    sys.exit(f"FAIL: типов тканей {len(set(tissues))}, нужно 5")

print("OK")
