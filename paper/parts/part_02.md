# Часть 2: с. 4–6

[с. 4]
www.nature.com/scientificreports/
RNA‑seq analysis. Trimming. We decided to apply the trimming procedure to eliminate adapter se-

quences present in our data and to improve read quality from the FASTQ files. Adapter removal and quality
trimming were carried out using Trimmomatic60 (v.0.35), Cutadapt61 (v.1.12) and BBDuk (v.Oct., 23, 2015), the
last included in the BBTools suite (https​://sourc​eforg​e.net/proje​cts/bbmap​/). In all cases only reads with a Phred
quality score > 20 and read length > 50 bp were selected for downstream analysis. The statistical calculations for
the comparisons among the trimming algorithms regarding the reads mapping rate and the surviving reads were
performed using the Kruskal–Wallis test, followed by the Dunn’s post-hoc test in the dunn.test62 (v.1.3.5) package
in ­R63 (v.3.5.1). Details about each trimming algorithm are given in Supplementary Note S1.
Alignment. We next evaluated the genome-based alignment methods Tophat216 (v.2.1.0), STAR​18 (v.2.5.3a)
and Hisat264 (v.2.0.0). The Human genome version GRCh37 (hg19) from ­Ensembl65 was used as the reference
genome. Transcriptome alignment methods represented by Bowtie266 (v.2.2.6) and STAR​were also tested against
the Ensembl (v82) transcriptome. Finally, hybrid methods, which combine both types of alignment, were represented by RUM19 (v.2.0.5_06) that used its own hg19 reference. A pre-processing step of adding Ns letters to the
incomplete reads was needed in RUM since paired reads with different lengths were not allowed. BAM files from
all the methods were sorted by read name and genome position using samtools67 (v.1.3.1) and unmapped reads
were discarded. The statistical calculations concerning the unmapped reads for each algorithm were made using
the Kruskal–Wallis test, followed by the Dunn’s post-hoc test in the dunn.test package. Details about each aligner
are given in Supplementary Note S1.
Counting and normalization. The number of alignments mapped to each gene was counted using as reference
a gene transfer format (GTF) file from Ensembl (v.82). We used six counting methods: Cufflinks68 (v.2.2.1),
eXpress69 (v.1.5.1), HTSeq70 (v.0.6.1p1) using both the Intersection-Strict and the Union approaches, RSEM27
(v.1.2.31) and Stringtie71 (v.1.3.3b). Gene expression values were represented using the normalization techniques
provided by each algorithm: Fragments per Kilobase of Mapped reads (FPKM), Transcripts per Million (TPM),
Trimmed Mean of M values (TMM from edgeR), Relative Log Expression (RLE from DESeq2), upper quartile
(UQ), coverage (Cov), estimated counts (Est_Counts) and effective counts (Eff_Counts). Gene level expression
values from RSEM and eXpress were obtained by summing the transcript-level estimates. Details about these
counting and normalization methods are given in Supplementary Note S1.
Pseudoalignment. These methods do not consider the classical alignment process and carry out alignment,
counting and normalization in one single step. We applied three commonly used pseudoaligners: Kallisto26
(v.0.43.1), Sailfish28 (v.0.9.2) and Salmon72 (v.0.8.2). Salmon was executed using the FMD and the quasi-mapping-based (QMB) indexing modes. Details about pseudoaligners are given in Supplementary Note S1.
Statistical approaches for pipeline precision analysis. Precision was tested on the six DMSO treated control samples (CLA-T0 and CLB-T0) considering each cell line separately (three replicates by cell line). For this purpose,
we considered the 107 HKg expressed in the six samples and in 32 healthy tissues previously reported (Supplementary Table S2). The median rank of the precision was calculated from the CoV of the 107 genes. The CoV
is defined as:

CoV =

MAD
Median

where MAD is the median absolute deviation for each gene, given by:

MAD = median(|Xi − median(X)|)
being X the expression value of the X gene and i the ith sample.
The choice of this parameter was justified because it is a non-parametric method that quantifies precision
independently of the measurement units, and it also considers the median variability through the MAD parameter, but in the same context as the median. Given this background, we first calculated the CoV value for every
single gene inside of each pipeline and each cell line. Secondly, we ranked CoV values in each of the 192 pipelines
generated from the combination of all the above algorithms and methods (Fig. 1). Thirdly, we computed the
median of the ranks for the 107 HKg in each pipeline and this value was considered as the precision index of
each pipeline. The lower the index, the more precise the pipeline was. This procedure is graphically described
in the upper panel of Fig. 2.
Statistical approaches for pipeline accuracy analysis. Accuracy was tested on the six DMSO treated samples
considering separately each cell line. We calculated the accuracy of each pipeline based on the 32 gene expression values from qRT-PCR, which is considered as a gold standard method due to its accuracy and ­sensitivity73.
The Pearson correlation coefficient (r) was used to assess the association between the RNA-seq data and the
qRT-PCR Ct values. Then, we calculated a rank for the correlation coefficient of each gene along each pipeline.
Finally, a median rank for each gene was calculated. The lower the rank, the more accurate the pipeline was. This
procedure is graphically described in the bottom panel of Fig. 2. All r correlation coefficients were calculated
using the corr.test function implemented in the psych74 R package (v.1.9.12.31).
TOP pipeline selection. The best pipelines were settled considering the summation of the precision and accuracy ranks inside each cell line, thus giving the same weight for both parameters. The lower the summation, the
Scientific Reports |
Vol:.(1234567890)

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

4

[с. 5]
www.nature.com/scientificreports/

Figure 3.  Experimental procedure. Two multiple myeloma cell lines (KMS12-BM [CLA] and JJN-3 [CLB]),
two drugs (Amiloride [T1] and TG003 [T2]), and dimethyl-sulfoxide (DMSO) (treatment 0 [T0]) were used
to conduct the RNA-seq and the qRT-PCR experiments. Control samples were used to carry out the raw gene
expression quantification study, whilst all the 18 samples were used to perform the differential gene expression
analysis.
better the performance of the pipeline was. All the comparisons involving the pipeline ranks were statistically
evaluated using the Kruskal–Wallis test followed by the Dunn’s post-hoc test. p-values were adjusted for multiple
testing by the Benjamini–Hochberg False Discovery Rate (FDR) ­procedure75. All these calculations were made
using the dunn.test package.
Differential gene expression methods for RNA‑seq. Seventeen DE methods were evaluated, as the result of the
combination of 11 DE algorithms and their different normalization options. Such methods can be divided into
three categories: (a) methods that assume a negative binomial distribution of data: baySeq76 (v.2.10.0), Cuffdiff77
(v.2.2.1), DESeq278 (v.1.16.1), EBseq79 (v.1.16.0) and the edgeR80 (v.3.18.1) generalized linear model (GLM), and
exact test variants; (b) methods that assume a log-normal distribution, like Ballgown81 (v.2.8.4) and the Trend and
Voom limma82 variants (v.3.32.10); and (c) non-parametric methods such as NOISeq48 (v.2.20.0) and SAMseq83
(available in the samr R package, v.2.0). Details about these methods are given in Supplementary Note S1. These
DE methods were executed using the top 10 pipelines regardless of the normalization step, because most of DE
methods only allow raw data as input. All these algorithms were assessed under five experimental contrasts, as
the result of the combination of the two cell lines with the two treatment options and DMSO (Fig. 3). Contrasts
were classified based on the number of DEGs detected. The sequence in descendent order of DEGs was: CLAT0 vs. CLB-T0 > CLA-T1 vs. CLA-T0 > CLA-T2 vs. CLA-T0 > CLB-T1 vs. CLB-T0 > CLB-T2 vs. CLB-T0. All the
above DEGs comparisons were carried out at 3 FDR cut-offs: FDR < 0.05, FDR < 0.01 and FDR < 0.001.
The similarity among the 17 methods was calculated in R through the dist function from the stats package
(v.3.5.1) using the Euclidean distance as distance measure and group average as linkage method. Dendrograms
were depicted using the dendextend package (v.1.13.3). The performance of the DE methods was determined
through the measurement of 7 diagnostic test parameters: sensitivity or true positive rate (TPR), specificity or
true negative rate (TNR), the positive predictive value (PPV), the negative predictive value (NPC), accuracy
(ACC), the area under the receiver operating characteristic curve (AUC) and the Matthews correlation coefficient (MCC). These parameters were calculated using as reference the qRT-PCR Benjamini-Hockberg’s adjusted
p-values in a two-sample t-test. The true positive, true negative, false positive and false negative estimates were
defined based on these qRT-PCR adjusted p-values and compared with the adjusted p-values obtained by each
RNA-seq DE method.

Equipment. All the computational procedures for raw gene expression quantification were performed in a
64-bit computer with 264 GB of RAM and 64 CPUs installed with a Linux system CentOS 6.9. The DGE quantification was carried out in a 64-bit computer with 64 GB of RAM and 24 CPUs with either a Linux system Ubuntu
14.04 or a Windows 10 system.
Scientific Reports |

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

5
Vol.:(0123456789)

[с. 6]
www.nature.com/scientificreports/

Ranking

Trimming method

Alignment method

Counting method

Normalization method

Median precision

Median accuracy

Median overall precision
and accuracy

1

Trimmomatic

RUM

HTSeq UNION

TMM

68.5

56

249

2

Trimmomatic

RUM

HTSeq INTER

TMM

68.5

57.75

252.5

3

BBDuk

RUM

HTSeq UNION

TMM

70

56.5

253

4

BBDuk

STAR​

HTSeq UNION

TMM

68

62

260

5

Cutadapt

TopHat2

HTSeq UNION

TMM

62.5

67.75

260.5

6

BBDuk

TopHat2

HTSeq UNION

TMM

63.5

68

263

7

BBDuk

HiSat2

HTSeq UNION

TMM

63.5

69.25

265.5

8

BBDuk

TopHat2

HTSeq INTER

TMM

62.5

70.5

266

9

Trimmomatic

STAR​

HTSeq UNION

TMM

69

64.75

267.5

10

Trimmomatic

STAR​

HTSeq INTER

TMM

63.5

71

269

Table 1.  Top 10 pipelines based on the overall precision and accuracy ranking. The lower the median overall
precision and accuracy the better the pipeline performance was.

Results

RNA‑seq workflow.

We tested the precision and accuracy of 192 RNA-seq pipelines in two independent and well-characterized MM cell lines at raw gene expression quantification level (RGEQ) (Fig. 3). These
192 pipelines are the result of the combination of different algorithms for trimming, alignment, counting and
normalization (Fig. 1, left panel). Comparisons were made using non-parametric statistics based on ranks and
median values as presented schematically in Fig. 2. We also tested 11 algorithms for DE combined with the normalization procedures available in their respective packages, which gave a total of 17 approaches (Fig. 1, right
panel).

Precision and accuracy of the evaluated pipelines. Precision was calculated using the 107 house-

keeping reference genes (HKg) (Supplementary Table S2) individually for each cell line, as described in “Methods” and in the upper panel of Fig. 2. This analysis showed that pipelines using pseudoalignment algorithms
were disproportionately represented among the top-ranked positions. Specifically, the most precise algorithm
was Salmon in conjunction with TPM normalization. The least precise pipelines were those in which the normalization process was less strict, such as the raw reads, Numreads, effective counts and estimated counts (Supplementary Table S4).
The accuracy of the 192 pipelines was estimated by testing by qRT-PCR 32 of the 107 genes used for the precision analysis (highlighted genes in Supplementary Table S2). Accuracy was calculated based on the qRT-PCR
results as described in the “Methods” section and in the bottom panel of Fig. 2. The top-ranked pipelines were
those in which HTSeq was the counting method and TMM was the normalization approach. It was also noticed
that the traditional alignment methods like RUM were situated in the top positions. The most accurate pipeline
used the Sailfish pseudoaligner. However, most pipelines that used a pseudoalignment method were ranked in
the lowest positions in the list for accuracy (Supplementary Table S5).
Finally, we generated a global ranking for all the tested pipelines, simultaneously considering the precision
and accuracy on RGEQ. The top 10 ranked pipelines are shown in Table 1, and the complete ordered list is
presented in Supplementary Table S6. Pipelines that used the counting algorithm HTSeq with its Union variant
and the TMM normalization were the most represented in the highest positions of the pipeline ranking. We
also observed that the RUM, Tophat2 and STAR​alignment algorithms appeared very frequently at the top positions of this list. We found no pattern for the trimming algorithms. These results clearly show that counting and
normalization methods are the most critical steps in the RNA-seq analysis process. Particularly, considering the
above results, we concluded that the combination of Trimmomatic + RUM + HTSeq Union + TMM was the most
precise and accurate pipeline.

Influence of individual RNA‑seq analysis steps on raw gene expression quantification. We
assessed the impact of each individual RNA-seq analysis step on the RGEQ. Thus, we performed this analysis at
five levels: trimming, alignment, counting and normalization, and the pseudoalignment.
Trimming. With respect to the effect of trimming, we explored the influence of three algorithms on RGEQ.
We found that the three algorithms displayed different patterns regarding the surviving paired reads and the
alignment rate. Thus, Cutadapt was the algorithm that produced the highest overall rate of surviving paired
reads (95.5%) and the lowest alignment rate (93.4%). On the other hand, BBDuk obtained the highest percentage of aligned reads (97.5%) (Supplementary Table S7). Despite these divergences, we did not find statistically
significant differences in RGEQ among the three algorithms (FDR > 0.05) (Fig. 4a and Supplementary Table S8).
Alignment. At this level, we evaluated five aligners that represent three alignment methodologies: to align
against the reference genome, to align against the reference transcriptome and a hybrid approach. This process
was evaluated on 156 pipelines because the 36 remaining pipelines were used to assess the pseudoalignment
performance. There were two parameters that could affected the RGEQ ranking at this level: the number of
Scientific Reports |
Vol:.(1234567890)

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

6

