---
type: Задача
status: выполнена
---
# Переписать tissues.py на pipeline

Запрос пользователя: «Перепиши src/tissues.py так, чтобы нормализация, отбор генов, фолды и обучение брались из src/pipeline.py (cpm_log, prepare_features, make_folds, fit_predict_proba), логику не копируй. checks/, METHOD.md, pipeline.py не меняй. После правки запусти python src/tissues.py и python3 src/check_result.py, покажи код возврата.»

Готово, когда: `python src/tissues.py` и `python3 src/check_result.py` завершаются без ошибок; `src/tissues.py` импортирует и использует `cpm_log`, `prepare_features`, `make_folds`, `fit_predict_proba` из `pipeline.py`; `checks/`, `METHOD.md`, `pipeline.py` не изменены.
Не входит: изменение `checks/`, `METHOD.md`, `pipeline.py`; изменение исходных данных.

## Итог
`src/tissues.py` переписан: нормализация (`cpm_log`), отбор генов (`prepare_features`), фолды (`make_folds`) и обучение (`fit_predict_proba`) импортируются из `src/pipeline.py` и используются без копирования логики. Локально остались только маркерные гены и генерация артефактов (barplot, radarplot, tissues_result.png/csv, method_list.md). `checks/`, `METHOD.md`, `pipeline.py` не изменены (git status — только `src/tissues.py` и результаты). Проверено: `python src/tissues.py` → accuracy 1.0000, код возврата 0; `python3 src/check_result.py` → OK, код возврата 0.

## Связи
- тип: [[Задача]]