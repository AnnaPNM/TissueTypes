# Часть 1: с. 1–3

[с. 1]
www.nature.com/scientificreports

OPEN

Systematic comparison
and assessment of RNA‑seq
procedures for gene expression
quantitative analysis
Luis A. Corchete1,2,3,4,5*, Elizabeta A. Rojas1,2,3, Diego Alonso‑López2, Javier De Las Rivas2,3,
Norma C. Gutiérrez1,2,3,5 & Francisco J. Burguillo4
RNA-seq is currently considered the most powerful, robust and adaptable technique for measuring
gene expression and transcription activation at genome-wide level. As the analysis of RNA-seq data
is complex, it has prompted a large amount of research on algorithms and methods. This has resulted
in a substantial increase in the number of options available at each step of the analysis. Consequently,
there is no clear consensus about the most appropriate algorithms and pipelines that should be used
to analyse RNA-seq data. In the present study, 192 pipelines using alternative methods were applied
to 18 samples from two human cell lines and the performance of the results was evaluated. Raw gene
expression signal was quantified by non-parametric statistics to measure precision and accuracy.
Differential gene expression performance was estimated by testing 17 differential expression
methods. The procedures were validated by qRT-PCR in the same samples. This study weighs up the
advantages and disadvantages of the tested algorithms and pipelines providing a comprehensive
guide to the different methods and procedures applied to the analysis of RNA-seq data, both for the
quantification of the raw expression signal and for the differential gene expression.
In recent years, RNA-seq has emerged as an alternative method to that of classic microarrays for transcriptome
­analysis1–4. Compared with microarrays, RNA-seq enables the study of novel transcripts and offers higher resolution, a better range of detection and lower technical ­variability5,6. RNA-seq also offers a high degree of agreement
with other techniques considered as the gold standard in transcriptomics such as qRT-PCR, both at absolute and
relative gene expression measurement ­levels7. All these facts have led to a great expansion of RNA-seq, becoming
the first choice in transcriptomic analysis for many scientists. Nevertheless, its widespread use has generated a
large amount of research on algorithms and methods that has eventually produced a lack of consensus about
how to analyse RNA-seq data. Thereby, in the last decade, many different algorithms and pipelines have been
proposed, but there is much debate about which approaches provide the most precise and accurate results. Thus,
further research to compare these methods remains necessary.
RNA-seq data analysis typically involves several steps: trimming, alignment, counting and normalization
of the sequenced reads, and, very often, differential expression (DE) analysis across conditions. Trimming is
used to increase reads mapping rate through the elimination of the adapter sequences and the removal of poorquality nucleotides. It must be employed non-aggressively, together with a wisely chosen read length, to avoid
unpredictable changes in gene ­expression8 and transcriptome ­assembly9. The beneficial effects of this process
have been evaluated in reference-based10,11 and de ­novo12,13 RNA-seq analyses. The alignment to a reference
genome or transcriptome is normally the second step in the RNA-seq workflow and has been widely evaluated
by many ­authors14–19. Once the reads have been mapped, they must be assigned to a gene or a transcript, in a
process known as counting or quantification. This is followed by a normalization procedure to remove possible
sequencing bias. Since counting followed by normalization is a crucial component of RNA-seq data analysis,
several methods have been developed and many comparative studies evaluating their suitability have been
­published20–28. The final step in most RNA-seq studies is DE analysis. The main concern at this point is how to
1

Hematology Department, University Hospital, 37007 Salamanca, Spain. 2Cancer Research Center (CiC‑IBMCC,
CSIC/USAL), Consejo Superior de Investigaciones Científicas (CSIC) and University of Salamanca (USAL),
37007 Salamanca, Spain. 3Institute of Biomedical Research of Salamanca (IBSAL), 37007 Salamanca,
Spain. 4Faculty of Pharmacy, University of Salamanca, 37007 Salamanca, Spain. 5Center for Biomedical Research
in Network of Cancer (CIBERONC), Salamanca, Spain. *email: lacorsan@usal.es
Scientific Reports |

(2020) 10:19737

| https://doi.org/10.1038/s41598-020-76881-x

1
Vol.:(0123456789)

[с. 2]
www.nature.com/scientificreports/

Figure 1.  RNA-seq analysis workflow. Left panel (1) represents the raw gene expression quantification
workflow. Every box contains the algorithms and methods used for the RNA-seq analysis at trimming,
alignment, counting, normalization and pseudoalignment levels. The right panel (2) represents the algorithms
used for the differential gene expression quantification. *HTSeq was performed in two modes: union and
intersection-strict. **EdgeR exact test, edgeR GLM and NOISeq have internally three normalization techniques
that were evaluated separately.
accurately detect differentially expressed genes (DEGs) between two or more conditions. For that purpose, a
large number of tools have been developed to facilitate the analysis, and thereby to construct lists of significant
­DEGs29–38. Due to the impact of this final step, the algorithms involved in DEG detection have been compared
in many publications where authors analyse the pros and cons of applying different ­algorithms23,29,39–48.
In this scenario, the major challenge in RNA-seq analysis is that the different steps must be sequentially
combined in a complete workflow or pipeline and users have to choose between many possible methodological
approaches and options. This shows the complexity of RNA-seq analysis and it is one of the most critical points
to obtain accurate results both at raw and DE expression levels. Many possible combinations of the current algorithms in RNA-seq have been comparatively analysed to help decide the best workflow, but their performance
remains under discussion. For example, Nookaew et al.49 investigated several workflows to pass from RNA-seq
reads to differential gene expression (DGE) in Saccharomyces cerevisiae. Seyednasrollah et al.43 compared software
packages for DE using two public data sets. Teng et al.50 suggested new metrics to assess differences between
competing pipelines and Williams et al.51 evaluated many workflows for DE in human monocytes sequenced as
51 bp single-end reads.
Facing with the problems described, the present work aims to evaluate 192 alternative methodological pipelines using Illumina HiSeq 2500 paired-end RNA-seq data measuring its accuracy and precision at raw gene
expression quantification level (RGEQ). These pipelines were constructed performing all the possible combinations of 3 trimming algorithms, 5 aligners, 6 counting methods, 3 pseudoaligners and 8 normalization
approaches. We also evaluated the performance of 17 DE methods using the results of the top 10 ranked RGEQ
pipelines. All these algorithms and methods were selected based on their popularity as they are used in dozens of
scientific publications and can be easily implemented into a pipeline (Fig. 1). As samples, we used data from two
multiple myeloma (MM) cell lines, treated with two different drugs as putative therapeutic molecules. A benchmark evaluation protocol for the accuracy and precision of each pipeline was set up, based in the experimental
testing of 32 genes by qRT-PCR and the detection of 107 house-keeping reference expressed genes (presented
in Fig. 2 and later described in more detail).

Methods

RNA‑seq data. We used two extensively characterized MM cell lines, KMS12-BM (cell line A [CLA])

and JJN-3 (cell line B [CLB]), as previously ­described52. Both cell lines were treated with Amiloride at 0.1 mM
(KMS12-BM) and 0.4 mM (JJN-3) for 24 h (treatment 1 [T1]), with TG003 at 50 µM (both KMS12-BM and
JJN-3) for 24 h (treatment 2 [T2]), and dimethyl-sulfoxide (DMSO) (treatment 0 [T0]) was used as negative
control. All the experiments were done by triplicate, making up a total of 18 samples. RNA was extracted using
Scientific Reports |
Vol:.(1234567890)

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

2

[с. 3]
www.nature.com/scientificreports/

Figure 2.  Benchmark procedure to evaluate precision and accuracy. Description of the procedure to evaluate
the precision (top) and the accuracy (bottom) in the RNA-seq analysis.
the RNeasy Plus Mini kit (QIAGEN, Hilden, Germany). The RNA integrity was assessed with an Agilent 2100
Bioanalyzer (Agilent Technologies). The corresponding 18 RNA libraries were constructed following the TruSeq
Strand-Specific RNA sequencing library protocol from Illumina (https​://suppo​rt.illum​ina.com/downl​oads/truse​
q_stran​ded_total​_rna_sampl​e_prepa​ratio​n_guide​_15031​048.html). The RNA libraries were sequenced in a
HiSeq 2500 system at Lifesequencing S.L. (Valencia, Spain). Paired-end reads of 101 base pairs (within a range
of 36,240,231–77,906,369 total reads), were generated using this platform (Supplementary Table S1). The quality
of the resulting sequences in FASTQ format was assessed using the FASTQC (v0.11.3) ­tool53.

Gene expression validation using qRT‑PCR. To determine qRT-PCR candidate genes we first selected

in our dataset 1181 genes expressed in 32 healthy tissues from the RNA-seq data published by Uhlen et al.54. We
then performed a filtering process removing all genes with < 4 gene expression units (alignments, FPKM, TPM,
etc. depending on the normalization method) in the six DMSO treated samples for each pipeline. This process
resulted in the choice of 107 genes that satisfied this filtering criterion in all the 192 pipelines (Supplementary
Table S2). As these 107 were constitutively expressed, both in 32 healthy tissues and in our two cell lines, we considered them as our gene expression reference set or housekeeping gene set (HKg). We selected from this first list
30 genes, 10 with the highest and 10 with the lowest median non-parametric coefficient of variation (CoV)55, and
10 random genes from the mid-area. Additionally, we selected two commonly used housekeeping genes in RNA
studies: GAPDH and ACTB, generating a second list of 32 genes. Starting from the same samples used in RNAseq, 1 µg of total RNA was reverse transcribed to cDNA using oligo dT from SuperScript First-Strand Synthesis
System for RT-PCR (Thermo Fisher Scientific). Taqman qRT-PCR mRNA assays (Applied Biosystems) were
carried out in duplicate on these 32 genes that are highlighted in the Supplementary Table S2. Oligonucleotide
probes used for the qRT-PCR assays are reported in the Supplementary Table S3.
To measure DGE by qRT-PCR we used the ΔCt method, calculated as:

Ct = CtControlgene −CtTargetgene
Three normalization approaches were tested: a) Endogenous control normalization, where endogenous control
was calculated as the mean of GAPDH and ACTB Ct values for each sample, b) Global median normalization, in
which the normalization factor was calculated using the median value for genes with Ct < 35 for each sample, and
c) Most stable gene, that was determined using the 4 algorithms available in the RefFinder webtool, BestKeeper56,
NormFinder57, Genorm58 and the comparative delta-Ct method59. The gene ECHS1 was ranked the top position
and it was considered as the most stable in our dataset. We detected bias on GAPDH and ACTB genes expression
due to the drug treatments. This bias was an under-expression of these genes in the treatment conditions as it is
shown in Supplementary Fig. S1. In view of this fact, we rejected the (a) normalization method. Regarding the
(b) and (c) methods, both performed equally well but the (c) method looked more robust than (b) as it captures
better the Ct values dispersion inside each sample. Finally, we chose the (b) method (Global median normalization) for our downstream analysis.

Scientific Reports |

(2020) 10:19737 |

https://doi.org/10.1038/s41598-020-76881-x

3
Vol.:(0123456789)

