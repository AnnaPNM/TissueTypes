# Применённые методы

## 1. Нормализация данных
RNA-seq каунты нормализованы в **log2(CPM+1)** (counts per million), чтобы
скомпенсировать разную глубину секвенирования образцов.

## 2. Отбор признаков
Выбраны **2000 генов** с наибольшей дисперсией log2(CPM) — они
несут основную тканеспецифичную информацию.

## 3. Классификация
- **RandomForest** (300 деревьев): accuracy = 0.9998 (5-fold CV)
- **Логистическая регрессия** (после StandardScaler): accuracy = 0.9998 (5-fold CV)

## 4. Определение типов тканей
Для каждого класса выделены топ-30 маркерных генов по разнице средних
log2(CPM) с остальными классами. По известным тканеспецифичным генам
(сверка с литературой) классы сопоставлены с тканями:

| Класс | Ткань | Ключевые маркеры |
|-------|-------|------------------|
| 1 | Skin | KRT14, KRT1, KRT5, KRT10, FLG, LORICRIN, DSP, DSG1 |
| 2 | Lung | SFTPB, SFTPA1, SFTPA2, SFTPC, NKX2-1, NAPSA, AGER |
| 3 | Brain | GFAP, PLP1, MOBP, SNAP25, SLC1A2, KIF5A, TUBB4A |
| 4 | Blood | HBB, HBA1, HBA2, HBD, S100A9, LYZ, MMP9 |
| 5 | Kidney | UMOD, AQP2, SLC12A1, SLC12A3, NPHS2, PAX2, PAX8, CDH16 |

## 5. Верификация
- Кросс-валидация (5-fold, стратифицированная) подтверждает точность ~100%.
- Маркерные гены каждого класса однозначно соответствуют известным
  тканеспецифичным генам.

## Источники
1. GTEx Consortium. The Genotype-Tissue Expression (GTEx) project.
   *Nature Genetics*, 2013.
2. Uhlen M. et al. Proteomics. Tissue-based map of the human proteome.
   *Science*, 2015.
3. Robinson M.D., Oshlack A. A scaling normalization method for differential
   expression analysis of RNA-seq data. *Genome Biology*, 2010.
4. Law C.W. et al. voom: precision weights unlock linear model analysis tools
   for RNA-seq read counts. *Genome Biology*, 2014.
5. Breiman L. Random Forests. *Machine Learning*, 2001.
