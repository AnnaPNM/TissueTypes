# Поиск по своей библиотеке: Tissue type identification via RNA-seq data

[источник] OpenAlex, запрос: <https://api.openalex.org/works?search=Tissue+type+identification+via+RNA-seq+data&per_page=50&filter=has_abstract%3Atrue&select=id%2Cdoi%2Cdisplay_name%2Cpublication_year%2Ccited_by_count%2Cabstract_inverted_index>, дата: 2026-10-21, работ: 20

## Какие есть пакеты для дифференциального анализа RNA-seq данных

| место | BM25                                                                       | вектор                                                                     |
|------------------------|------------------------|------------------------|
| 1     | Comparative analysis of RNA-Seq alignment algorithms and the (W2096173343) | A survey of best practices for RNA-seq data analysis (W2236822143)         |
| 2     | Statistical Modeling of RNA-Seq Data (W2068568257)                         | Single-Cell RNA-Seq Technologies and Related Computational D (W2926010705) |
| 3     | VIPER: Visualization Pipeline for RNA-seq, a Snakemake workf (W2799170276) | Detecting differential usage of exons from RNA-seq data (W2119204091)      |

2 и 3 место BM25 - по теме, 1 под вопросом, потому что там упор не на диф анализ, а на выравнивание

1 место у вектора - самое релевантное из всех источников\
Но при этом 2 и 3 немного отходят от темы

## Какие есть открытые базы данных с маркерными генами для разных организмов

| место | BM25 | вектор                                                                     |
|------------------------|------------------------|------------------------|
| 1     | ---  | Verifying explainability of a deep learning tissue classifie (W3125482371) |
| 2     | ---  | Gene Mapping via Bulked Segregant RNA-Seq (BSR-Seq) (W2072715465)          |
| 3     | ---  | EPISCORE: cell type deconvolution of bulk tissue DNA methylo (W3083730060) |

BM25 ничего не нашел (?)

Векторный поиск нашел статьи, в которых используются открытые базы данных

## Какие есть методы для определения маркерных генов

| место | BM25 | вектор                                                                     |
|------------------------|------------------------|------------------------|
| 1     | ---  | EPISCORE: cell type deconvolution of bulk tissue DNA methylo (W3083730060) |
| 2     | ---  | A survey of best practices for RNA-seq data analysis (W2236822143)         |
| 3     | ---  | Single-Cell RNA-Seq Technologies and Related Computational D (W2926010705) |

BM25 снова ничего не нашел

Векторный поиск нашел статью по деконволюции - это буквально процесс определения тканей по экспрессии генов (как в запросе)

## Может ли набор генов быть маркерным для нескольких тканей и как отличить эти ткани в таком случае

| место | BM25 | вектор                                                                     |
|------------------------|------------------------|------------------------|
| 1     | ---  | RNA-seq: from technology to biology (W2004419591)                          |
| 2     | ---  | Verifying explainability of a deep learning tissue classifie (W3125482371) |
| 3     | ---  | Identification of novel transcripts in annotated genomes usi (W2098832713) |

BM25 снова ничего не нашел

Вектор выдал то, что может помочь с ответом на вопрос (статьи релевантны)

## Где применяется определение типа ткани по RNA-seq данным

| место | BM25                                                                       | вектор                                                                |
|------------------------|------------------------|------------------------|
| 1     | Comparative analysis of RNA-Seq alignment algorithms and the (W2096173343) | RNA-seq: from technology to biology (W2004419591)                     |
| 2     | Statistical Modeling of RNA-Seq Data (W2068568257)                         | Detecting differential usage of exons from RNA-seq data (W2119204091) |
| 3     | VIPER: Visualization Pipeline for RNA-seq, a Snakemake workf (W2799170276) | A survey of best practices for RNA-seq data analysis (W2236822143)    |

BM25 нашел больше алгоритмов, а векторный поиск - статьи/обзоры непосредственно с применением описанного процесса (как в запросе)
