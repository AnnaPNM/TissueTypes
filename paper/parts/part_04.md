# Часть 4: с. 11–12

[с. 11]
www.nature.com/scientificreports/

Figure 7.  Summary of the performance of the RNA-seq gene differential expression analysis methods. This
graph includes three experimental approaches for the 17 methods: performance by number of DEG scenario,
performance by statistical significance cut-off and overall performance.
All these performance analyses are summarised jointly considering all the previous approaches in Fig. 7. We
highlight that the most balanced method was limma trend since, unlike the other 16 methods, it was ranked
among the 8 best methods in all the approaches. Other methods with good performance in all the scenarios were
baySeq, NOISeq FPKM, limma voom and some variants of edgeR GLM.
To complete the DE analysis, we estimated the influence of normalization on DEGs detection by the DE
methods that admitted multiple normalization approaches. We observed subtle differences in DEGs detection
between the normalization approaches with regard to AUC and ACC, particularly in some analysis scenarios of
the NOISeq algorithm, in which the FPKM normalization notably surpassed the other methods (Supplementary
Fig. S7).

Discussion

The findings of our study propose an optimal workflow for RNA-seq gene expression data analysis, based on the
performance of multiple algorithms and pipelines freely available to the scientific community.
On a first stage, we performed a comparison of 192 pipelines to quantify raw gene expression. According to
our results, the pipeline that obtained the best precision and accuracy rankings, and therefore it was the best
combination of algorithms, was Trimmomatic + RUM + HTSeq Union + TMM. Despite this fact, the combination of methods that made the difference were HTSeq at the counting level and TMM at normalization level, as
they were part of 10 out of the top 10 pipelines. On the other hand, the pipelines consistent on raw alignments,
Stringtie’s coverage and Salmon’s NumReads occupied most of the last positions in the ranking, confirming the
importance of a proper combination of counting and normalization methods.
In this work, we also compared the most commonly used methods over the stages of RNA-seq analysis: trimming, alignment, counting and normalization. Interestingly, we did not find substantial differences in RGEQ
related to the trimming algorithm or the aligner used in the RNA-seq data analysis. This last finding is consistent
with the recent work of Schaarschmidt et al. in Arabidopsis thaliana ­samples84, showing that the seven analysed
aligners could be equally applied in RNA-seq data analysis. On the other hand, we found that the counting and
normalization step was shown to be a critical step in our RNA-seq data analysis. This has also been reported by
Robert et al.85, who demonstrated that the counting step may overestimate or underestimate the level of gene
expression. We also detected a preponderance of the HTSeq counting method (under the Union and Intersectionstrict variants) in the top-ranked positions in our pipeline list. With respect to data normalization, pipelines based
on the TMM method outperformed those using RLE, FPKM, TPM and all the other normalization approaches
tested. This finding is consistent with those of Wu et al.21 and Maza et al.20, but is not in agreement with those
of Li et al.22, who concluded that gene normalization methods did not improve the gene expression calculations
provided by the raw counts.
We also measured the performance of the pseudoalignment on RGEQ. Recently, pseudoaligners have entered
the scene as an alternative to the classic alignment algorithms. One of their main advantages is that they can
carry out the processes of alignment, counting and normalization in a single step. We have seen that pipelines
Scientific Reports |

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

11
Vol.:(0123456789)

[с. 12]
www.nature.com/scientificreports/
based on pseudoaligners have a good precision in gene expression estimation, but their accuracy is inferior to
that of the classic aligners. In any case, pseudoaligners could be a good alternative as RNA-seq data exploration
tools since their execution time is faster than the conventional alignment methods.
Regarding RNA-seq DGE quantification, we established that the performance of the methods depended on
the number of DEGs present under the contrasted experimental conditions and the FDR cut-off employed to
determine statistical significance. Our systematic analysis revealed that limma trend obtained the best results in
terms of performance, closely followed by limma voom, NOISeq FPKM, baySeq, and some derivations of edgeR.
It is interesting that even though the assumptions about the underlying distribution differ, the performance of
these methods is comparable, as was also described by Rapaport et al.39 with other sets of methods. Of note, in
spite of their similar performance, we detected differences in behavioural patterns for these methods. Whilst
limma trend, limma voom and edgeR showed a good balance among the seven diagnostic test parameters analysed, baySeq and NOISeq tended to be biased in favour of some of them. In the case of baySeq, we found that its
performance was sustained by the excellent specificity of this method, however, we detected a substantial lack
of sensitivity. With respect to NOISeq, its performance was supported by its good sensitivity and precision, but
we found low levels of specificity. A similar pattern of unequal performance between diagnostic test parameters
for NOISeq was also reported in the work of Williams et al.51, who concluded that NOISeq exhibited the lowest
recall (sensitivity), but agreed with our work in the high precision of this algorithm.
Another important aspect of the DE analysis is the possible influence of normalization on the results. Normalization is a critical step in RNA-seq DGE analysis since it enables expression levels to be compared. According
to our results, the choice of the normalization method had little effect on DEGs detection compared to RGEQ.
These findings are supported by the research of Assefa et al.86 and Seyednasrollah et al.43, who reported a considerable overlap among DE methods regardless of the normalization approach used. However, our results are in
contrast with those of Dillies et al.25, Bullard et al.87 and Zyprych-Walczak et al.88, possibly because these authors
evaluated several normalization approaches regardless of the DE method, while in our study we only assessed
the normalization provided by each DE method. It is also of note that Ballgown and SAMseq performed worse
than others when the compared conditions had a low number of DEGs. It must be borne in mind that this bad
behaviour may lead to a high frequency of false-negative deregulated genes in the data analysis.

Data availability

RNA-seq data used in this study are deposited in the Gene Expression Omnibus (GEO, https​://www.ncbi.nlm.
nih.gov/geo/) under the accession number GSE95077.
Received: 15 May 2020; Accepted: 3 November 2020

References

1. Garber, M., Grabherr, M. G., Guttman, M. & Trapnell, C. Computational methods for transcriptome annotation and quantification
using RNA-seq. Nat. Methods 8, 469–477 (2011).
2. Xuan, J., Yu, Y., Qing, T., Guo, L. & Shi, L. Next-generation sequencing in the clinic: promises and challenges. Cancer Lett. 340,
284–295 (2013).
3. Finotello, F. & Di Camillo, B. Measuring differential gene expression with RNA-seq: challenges and strategies for data analysis.
Brief Funct. Genomics 14, 130–142 (2015).
4. Han, Y., Gao, S., Muegge, K., Zhang, W. & Zhou, B. Advanced applications of RNA sequencing and challenges. Bioinform Biol.
Insights 9, 29–46 (2015).
5. Perkins, J. R. et al. A comparison of RNA-seq and exon arrays for whole genome transcription profiling of the L5 spinal nerve
transection model of neuropathic pain in the rat. Mol. Pain 10, 7 (2014).
6. Zhao, S., Fung-Leung, W. P., Bittner, A., Ngo, K. & Liu, X. Comparison of RNA-Seq and microarray in transcriptome profiling of
activated T cells. PLoS ONE 9, e78644 (2014).
7. SEQC/MAQC-III Consortium. A comprehensive assessment of RNA-seq accuracy, reproducibility and information content by
the Sequencing Quality Control Consortium. Nat. Biotechnol. 32, 903–914 (2014).
8. Williams, C. R., Baccarella, A., Parrish, J. Z. & Kim, C. C. Trimming of sequence reads alters RNA-Seq gene expression estimates.
BMC Bioinformatics 17, 103 (2016).
9. Macmanes, M. D. On the optimal trimming of high-throughput mRNA sequence data. Front. Genet. 5, 13 (2014).
10. Chen, C., Khaleel, S. S., Huang, H. & Wu, C. H. Software for pre-processing Illumina next-generation sequencing short read
sequences. Source Code Biol. Med. 9, 8–0473–9–8. eCollection 2014 (2014).
11. Del Fabbro, C., Scalabrin, S., Morgante, M. & Giorgi, F. M. An extensive evaluation of read trimming effects on Illumina NGS data
analysis. PLoS ONE 8, e85024 (2013).
12. Garg, R., Patel, R. K., Tyagi, A. K. & Jain, M. D. novo assembly of chickpea transcriptome using short reads for gene discovery and
marker identification. DNA Res. 18, 53–63 (2011).
13. Mbandi, S. K., Hesse, U., Rees, D. J. & Christoffels, A. A glance at quality score: implication for de novo transcriptome reconstruction of Illumina reads. Front. Genet. 5, 17 (2014).
14. Borozan, I., Watt, S. N. & Ferretti, V. Evaluation of alignment algorithms for discovery and identification of pathogens using RNASeq. PLoS ONE 8, e76935 (2013).
15. Yang, C., Wu, P. Y., Tong, L., Phan, J. H. & Wang, M. D. The impact of RNA-seq aligners on gene expression estimation. ACM BCB
2015, 462–471 (2015).
16. Kim, D. et al. TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome
Biol. 14, R36 (2013).
17. Yang, C., Wu, P. Y., Phan, J. H. & Wang, M. D. The impact of RNA-seq alignment pipeline on detection of differentially expressed
genes. IEEE Glob. Conf. Signal. Inf. Process. 2012, 1376–1379 (2014).
18. Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21 (2013).
19. Grant, G. R. et al. Comparative analysis of RNA-Seq alignment algorithms and the RNA-Seq unified mapper (RUM). Bioinformatics
27, 2518–2528 (2011).
20. Maza, E., Frasse, P., Senin, P., Bouzayen, M. & Zouine, M. Comparison of normalization methods for differential gene expression
analysis in RNA-Seq experiments: a matter of relative size of studied transcriptomes. Commun. Integr. Biol. 6, e25849 (2013).

Scientific Reports |
Vol:.(1234567890)

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

12

