# Часть 3: с. 7–10

[с. 7]
www.nature.com/scientificreports/

Figure 4.  Influence of the algorithms on the RNA-seq raw gene expression quantification. Box-plot analysis of the 192
pipelines grouped by the algorithms used at each step of the procedure: (a) trimming algorithms, (b) alignment algorithms (c)
counting methods, (d) normalization methods and (e) pseudoalignment algorithms. Coloured boxplots represent the scaled
values (between 1 and 100) of the summation of the precision and accuracy ranking of the 192 pipelines before (blue) and
after (green) the removal of pipelines that used raw reads, effective counts, estimated counts and coverage, which produced
a bimodal data distribution. The red diamond represents the mean of the ranking reached by the pipelines that use the
respective method or algorithm. The asterisks indicate the significance of the post-hoc Dunn’s test: *p < 0.05, **p < 0.01 and
***p < 0.001. Comparisons without asterisk are statistically insignificant (p > 0.05). Asterisks in (d) correspond to the p-values
of the most significant method (TMM) against the other methods.
Scientific Reports |

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

7
Vol.:(0123456789)

[с. 8]
www.nature.com/scientificreports/
concordant alignments and the unmapped reads output by each aligner. We found huge differences for these
two parameters among the five aligners. STAR​ was the algorithm with a higher number of concordant alignments (median = 93%, MAD = 1.0), closely followed by TopHat2 (Median = 90.1%, MAD = 2.0). On the other
hand, Bowtie2 only reached a 41.5% (MAD = 1.7) of concordant alignments, showing statistically significant
differences with the other four aligners (FDR < 0.05) (Supplementary Fig. S2). Regarding to the unmapped reads,
HiSat2, and STAR​outperformed the other aligners (FDR < 0.05), both producing less than 5% of unmapped
reads, although HiSat2 was particularly effective with only 1.6% (MAD = 0.3) unmapped reads (Supplementary
Fig. S2).
In view of this, we tested the influence of these aligners on the RGEQ and, despite their particularities, we did
not find any statistically significant differences (FDR > 0.05) between the median ranks of the pipelines involving
these alignment algorithms (Fig. 4b and Supplementary Table S9). However, we found a bimodal distribution in
the RGEQ ranking values for these algorithms. We then investigated the cause of this bimodal distribution and
we discovered that it was associated with the method employed for normalization. Therefore, the pipelines that
used raw reads, effective counts, estimated counts and coverage, reached the lowest ranks. A reanalysis removing
pipelines with these normalization approaches revealed that Bowtie2 pipelines reached poorer ranks of RGEQ
than pipelines that employed other alignment methods (Fig. 4b).
Counting and normalization. After the alignment step, we evaluated the influence of the counting methods
on RGEQ. Pipelines based on Cufflinks and RSEM reached the highest positions in the pipeline ranking followed by HTSeq and Stringtie based pipelines. These last two methods also obtained high-ranking positions, but
they showed a bimodal distribution in their rank values that was eventually translated in a greater variability,
similarly to that observed in the alignment step (Fig. 4c). The removal of the pipelines that used raw reads,
effective counts, estimated counts and coverage as normalization approach, caused the clearance of the bimodal
distributions and revealed a better performance of HTSeq pipelines (HTSeq-Inter and HTSeq-Union, Fig. 4c). It
also should be noted that the eXpress pipelines showed the poorest rates for the sum of precision and accuracy
(Fig. 4c and Supplementary Table S10).
Next, we evaluated the influence of normalization methods after the counting procedure. Pipelines using the
Trimmed Mean of M values (TMM) method performed the best on RGEQ. Other normalization methods, such as
Relative Log Expression (RLE) that was second best, and Transcript Per Million (TPM) or Fragments Per Kilobase
of Mapped reads (FPKM), also reached high-ranking positions for precision and accuracy in RGEQ analysis
(Fig. 4d and Supplementary Table S11). Considering these results, the normalization procedure proved to be
an essential step in RGEQ, since we detected the higher statistically significant differences among the pipeline
ranks based on the normalization method used.
Quantification by pseudoalignment. The RGEQ analysis was completed with the evaluation of three pseudoalignment algorithms considering their different alignment modes. These algorithms have a major advantage over
traditional alignment methods as they carry out alignment and counting in a single step. Regarding the performance of these algorithms on RGEQ, we found no statistically significant differences (FDR > 0.05) between their
cumulative ranks of precision and accuracy (Fig. 4e and Supplementary Table S12). When we compared these
pseudoalignment algorithms with traditional aligners, we found a similar performance on RGEQ (FDR > 0.05)
among the algorithms tested from both methodologies (Supplementary Table S13).

Performance of differential gene expression methods. We tested 17 DE methods obtained from the

combination of the different DE and normalization approaches. They were tested under six experimental conditions with three biological replicates per condition. Comparisons under these conditions are explained in Fig. 3.
The efficiency of the DE procedures was checked at three commonly used FDR cut-offs: FDR < 0.05, FDR < 0.01
and FDR < 0.001.
First, we compared the detection ability among the DE methods. There was great homogeneity among the
methods tested in the comparisons with a large number of gene expression changes. However, this homogeneity
decreased as the compared experimental conditions became more similar (Fig. 5). Consistent with the findings of
Seyednasrollah et al.43, we discovered that Cuffdiff was generally the algorithm that detected the smallest number
of gene expression changes, whilst SAMseq was the method that detected the largest number of changes in all
compared conditions. It was also observed that SAMseq itself and Ballgown lost detection power at FDR < 0.01
and FDR < 0.001 levels. BaySeq, EBSeq and Cuffdiff also lost detection power when the contrasted conditions were
more similar. On the other hand, methods derived from edgeR, limma and NOISeq were stable in all comparisons
at the three FDR levels analysed (Fig. 5a–e). The similarity analysis through the Euclidean distance among these
17 methods detected a high correspondence among edgeR, limma and DESeq2 in all scenarios. However, methods
such as SAMseq, Cuffdiff, EBSeq and Ballgown showed greater distances than most of the other methods. NOISeq,
meanwhile, was the most variable method regarding the DEG scenario (Supplementary Fig. S3).
We also analysed 7 diagnostic test parameters for each DE method using qRT-PCR as described in the “Methods” section. The individual analysis of each parameter revealed substantial differences among the 17 methods
(Fig. 6). Particularly, regarding the MCC (Fig. 6a), which evaluates the quality of the classification achieved
by each method, some DE methods such as NOISeq, Cuffdiff, baySeq and some of the edgeR variants showed a
moderate and positive relationship (MCC > 0.3) between their results and the true model. Concerning ACC and
AUC (Fig. 6b,c), the most accurate methods, that is, those with a higher ACC, were limma trend, limma voom
and baySeq at FDR < 0.001 (ACC = 0.78), while the method with a greater AUC, or in other words, the method
with the highest discrimination capacity, was baySeq at FDR < 0.001 (AUC = 0.81).

Scientific Reports |
Vol:.(1234567890)

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

8

[с. 9]
www.nature.com/scientificreports/

Figure 5.  Differential expression detection. Number of differentially expressed genes (DEGs) detected by the
17 methods of differential expression at three FDR cut-offs: 0.05, 0.01 and 0.001. Panels represent different
group comparisons in descending order based on the number of DEGs. (a) KMS12-BM (CLA) + DMSO (T0)
vs. JJN-3 (CLB) + DMSO (T0) (b) KMS12-BM (CLA) + Amiloride (T1) vs. KMS12-BM (CLA) + DMSO (T0)
(c) KMS12-BM (CLA) + TG003 (T2) vs. KMS12-BM (CLA) + DMSO (T0) (d) JJN-3 (CLB) + Amiloride (T1) vs.
JJN-3 (CLB) + DMSO (T0) (e) JJN-3 (CLB) + TG003 (T2) vs. JJN-3 (CLB) + DMSO (T0).
Next, we carried out a performance analysis considering simultaneously 7 parameters in three experimental
approaches: performance by number of DEG scenario, performance by statistical significance cut-off and overall
performance. Regarding the performance by the number of DEGs, we observed great variability among the 17
DE methods, so that the first ranking position was reached by different methods in the five analysed scenarios.
In this way, methods such as EBseq, DESeq2 or baySeq showed a highly variable behaviour depending on the
number of DEGs of the analysed scenario (Supplementary Fig. S4). On the other hand, the analysis of statistical
significance revealed some methods with fairly high and stable performance ranks: limma trend, limma voom
and baySeq (with also good ranks observed for edgeR exact and NOISeq, but suspected of greater variability
depending on whether they used UQ, FPKM or TMM) (Supplementary Fig. S5).
To carry out the third approach, we considered both the five DEG scenarios and the three statistical significance levels. Taken together, baySeq and NOISeq UQ were the methods that showed a better behaviour even
though both classified between the 9th and the 17th place in at least four of the approaches tested (Supplementary
Fig. S6). These two methods were closely followed by limma trend, limma voom, and edgeR GLM.

Scientific Reports |

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

9
Vol.:(0123456789)

[с. 10]
www.nature.com/scientificreports/

Figure 6.  Analysis of performance of the 17 differential gene expression methods through the measurement of
7 diagnostic test parameters. (a) Matthews correlation coefficient (MCC), (b) accuracy (ACC), (c) area under
the ROC curve (AUC), (d) positive predictive value (PPV), (e) negative predictive value, (f) true positive rate
(TPR), and (g) true negative rate (TNR). Performance was measured at three FDR cut-off levels: FDR < 0.05,
FDR < 0.01 and FDR < 0.001 for the 17 methods.

Scientific Reports |
Vol:.(1234567890)

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

10

