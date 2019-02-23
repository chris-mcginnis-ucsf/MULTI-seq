# MULTI-seq Sample Classification
MULTI-seq is a method for single-cell RNA sequencing sample multiplexing (for more information, check out our preprint: https://www.biorxiv.org/content/early/2018/08/08/387241). MULTI-seq is methodologically analogous to the Cell Hashing (Stoeckius et al., 2018, Genome Biology) and Click-Tags (Gehring et al., 2018, bioRxiv) except we utilize lipid- and cholesterol-modified oligonucleotides to label live-cell and nuclear membranes.

Here, we provide source code for two distinct stages of MULTI-seq data processing:

1. FASTQ conversion to a MULTI-seq sample barcode UMI count matrix. Think of this pipeline as 'CellRanger' for sample barcode data, as it (1) Splits raw FASTQs into cell barcode, UMI, and sample barcode sequences, (2) Removes reads that do not align with >1 mismatch to any MULTI-seq sample barcode reference sequence, (3) Removes reads representing duplicated UMIs for each cell, and (4) Convert this parsed read table to a sample barcode UMI count matrix. This count matrix can be used as the input for the MULTI-seq sample classification workflow (discussed below), or alternative classification strategies (Seurat, DemuxEM, etc.)

![alternativetext](/Figures/MULTIseq_Alignment_2.png)

2a. Sample classification. The MULTI-seq sample classification workflow expands on concepts borrowed from Cell-Hashing (Stoeckius et al., 2018) and Perturb-seq (Adamson et al., 2016; Dixit et al., 2016). As in Perturb-Seq, we model the probability density function for each sample barcode UMI distribution using Guassian-kernel density estimation. We then define local maxima corresponding to positive cells for each sample, as well as background cells. We then perform a inter-maxima quantile sweep to find the threshold for each barcode that results in the maximum number of singlet classifications. Using these thresholds, we then classify cells according to the number of barcode thresholds it surpasses (as in Cell Hashing) -- i.e., 0 thresholds = Negative, 1 threshold = Singlet, >1 threshold = Doublet/Multiplet.

![alternativetext](/Figures/MULTIseq_ClassificationWorkflow.png)

2b. We then utilize these initial classification results to 'rescue' a subset of previously-unclassified cells using the strategy employed in Cell Hashing. Specifically, we compute the classification stability (CS) for all negative cells. CS is the number of quantiles over which a given cell surpassed a single threshold. The thought is that putatively-reclassifiable cells will differ from true negatives and doublets in this way, as putative singlets will have higher CS than true negatives and doublets. We then use existing classification results as 'ground-truth' for initializing cluster centers during k-means clustering. We then compute the matching rate between k-means results and 'ground-truth' classifications, and compare this matching rate for negative cells binned by CS values. The CS value at which negative and 'ground-truth' matching rates deviate is the point at which negative cells can no longer be confidently reclassified.

![alternativetext](/Figures/MULTIseq_NegativeCellReclassification.png)

# Installation
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

# Tutorial: 96-plex HMEC sample multiplexed scRNA-seq

# Referencens
1. Stoeckius M, Zheng S, Houck-Loomis B, Hao S, Yeung BZ, Smibert P, Satija R. Cell "hashing" with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. 2017. Preprint. bioRxiv doi: 10.1101/237693.
2. Adamson B, Norman TM, Jost M, Cho MY, Nuñez JK, Chen Y, et al. A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response. Cell. 2016; 167(7):1867-82.e21.
3. Dixit A, Parnas O, Li B, Chen J, Fulco CP, Jerby-Arnon L, et al. Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens. Cell. 2016; 167(7):1853-66.e17.
4. Gaublomme JT, Li B, McCabe C, Knecht A, Drokhlyansky E, Van Wittenberghe N, Waldman J. Nuclei multiplexing with barcoded antibodies for single-nucleus genomics. 2018. Preprint. bioRxiv doi: 10.1101/476036.
