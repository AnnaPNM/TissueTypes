# STATUS.md — где работа сейчас

## Где мы
Все 5 зашифрованных типов расшифрованы:
- Тип 1 = **Skin** (принят пользователем ранее)
- Тип 2 = **Lung** (лёгкое)
- Тип 3 = **Brain** (мозг)
- Тип 4 = **Blood** (кровь)
- Тип 5 = **Kidney** (почка)

`src/tissues.py` переписан на использование `src/pipeline.py`
(`cpm_log`, `prepare_features`, `make_folds`, `fit_predict_proba`); логика не копируется.
`python src/tissues.py` → accuracy 1.0000, EXIT=0; `python3 src/check_result.py` → OK, EXIT=0.
`checks/`, `METHOD.md`, `pipeline.py` не изменены.

## Дальше
Задача «Переписать tissues.py на pipeline» выполнена — ожидает принятия пользователем.
Задача «Определить ткани для типов 2, 3, 4, 5» выполнена — ожидает принятия.
Задача «Собрать литературу по теме» — открыта.

## Действует
- Исходные данные в `data/` не изменять.
- `results/tissues_result.csv` создан скриптом `src/tissues.py`.
- Граф проекта: заметки классов в `ontology/`, заметки объектов в `notes/`.