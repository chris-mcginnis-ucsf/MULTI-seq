~~ Announcement (12/06/21) ~~
MULTI-seq reagents are now available for purchase from Millipore-Sigma: https://www.sigmaaldrich.com/US/en/product/sigma/lmo001

~~ Announcement (11/24/21) ~~
I'm now a postdoc at Stanford and my UCSF email will be decommissioned soon. I also only check my github repos about once per month, so please reach out directly at cmcginni@stanford[dot]edu if you run into any issues. 

# deMULTIplex (code written by Chris McGinnis)
deMULTIplex is an R package containing the companion software for our method for single-cell RNA sequencing sample Multiplexing Using Lipid-Tagged Indices: MULTI-seq (for more information, check out our Nature Methods manuscript: https://www.ncbi.nlm.nih.gov/pubmed/31209384)

MULTI-seq is methodologically analogous to the Cell Hashing (Stoeckius et al., 2018, Genome Biology) and Click-Tags (Gehring et al., 2019, NBT) except we utilize lipid- and cholesterol-modified oligonucleotides to rapidly and non-perturbatively label live-cell and nuclear membranes.

![alternativetext](/Figures/MULTIseq.overview.png)

We have signed an exclusive license with Millipore-Sigma for selling MULTI-seq reagents -- reagents will be available for direct purchase by the beginning of 2022, but you can get early access by reaching out to James Hoberg (james.hoberg@milliporesigma[dot]com). We are still providing LMOs through formal collaborations -- if this route is of interest, please reach out to Zev Gartner directly (zev.gartner@ucsf[dor]com).

Download example dataset for trying deMULTIplex here: https://stanfordmedicine.box.com/v/deMULTIplexdata

NOTE: This 8-donor PBMC experiment is different from the HMEC data used in the tutorial. I opted to upload these data because (i) the MULTI-seq barcoding work very well and there were no missing barcodes, (ii) negative-cell reclassification was unncessary because >95% of the cells were classified, and (iii) the files are much smaller. 

## Updates
(07/01/2020): Added link to download example dataset -- 8-donor PBMC experiment.

Version 1.0.2 (05/09/2019):
* Replaced 'chemistry' argument in 'MULTIseq.preProcess' and 'MULTIseq.preProcess_allCells' with 'cell', 'umi', and 'tag' arguments for more flexible read parsing (per the suggestion from @omansn). Each argument takes a length-2 numerical vector specifying the beginning and ending positions of the cell barcode (e.g., cell=c(1,16) for 10X), UMI (e.g., umi=c(17,28) for 10X V3 chemistry), and sample tag (e.g., tag=c(1,8) for MULTI-seq). Note: CellID and UMI sequences are assumed to be on R1, while sample tag sequences are assumed to be on R2.

Version 1.0.1 (03/13/2019):
* Added 'chemistry' arguments to 'MULTIseq.preProcess' and 'MULTIseq.preProcess_allCells' for 10X V2 or V3 chemistry. Argument changes the UMI position specification (e.g., R1 17-26 for V2, R1 17-28 for V3).

## Installation (in R/Rstudio)
devtools::install_github('chris-mcginnis-ucsf/MULTI-seq')

## Dependencies
deMULTIplex requires the following R packages: 
* KernSmooth (2.23-15) 
* reshape2 (1.4.3) 
* rTSNE (0.15) 
* stringdist (0.9.5.1)
* ShortRead (1.40.0)
* NOTE: These package versions were used in the bioRxiv paper, but other versions may work, as well.

# deMULTIplex Overview
## MULTI-seq sample barcode FASTQ alignment
R functions found in 'MULTIseq.Align.Suite.R' can be used to convert MULTI-seq sample barcode FASTQ files into a sample barcode UMI count matrix. Think of this pipeline as 'CellRanger' for sample barcode data. This pre-processing pipeline can be split into four distinct steps:
1. Split raw FASTQs into cell barcode, UMI, and sample barcode sequences associated with a user-defined set of cell barcodes
2. Remove reads that do not align with >1 mismatch to any MULTI-seq sample barcode reference sequence, 
3. Remove reads representing duplicated UMIs on a cell-by-cell basis
4. Convert this parsed read table to a sample barcode UMI count matrix. 

![alternativetext](/Figures/MULTIseq_Alignment.png)

This sample barcode UMI count matrix can be used as the input for the MULTI-seq sample classification workflow (discussed below), or alternative classification strategies (Seurat, DemuxEM, etc.).

## MULTI-seq sample classification workflow
R functions found in 'MULTIseq.Classification.Suite.R' can applied to MULTI-seq sample barcode UMI count matrices to classify cells into sample groups. The MULTI-seq sample classification workflow builds upon concepts borrowed from Cell-Hashing (Stoeckius et al., 2018) and Perturb-seq (Adamson et al., 2016; Dixit et al., 2016) and can be split into five distinct steps:
1. Model the probability density function for each normalized sample barcode UMI distribution using Gaussian-kernal density estimation (as in Perturb-Seq)
2. Define local maxima corresponding to positive cells (highest maxima) and background cells (mode)
3. Define barcode-specific thresholds across an inter-maxima quantile sweep (0.2-0.99), compute the proportion of classified singlets for each quantile
4. Using barcode-specific thresholds generated by the quantile that maximizes the proportion of singlets, classify cells according to the number of barcode thresholds it surpasses (as in Cell Hashing) -- i.e., 0 thresholds = Negative, 1 threshold = Singlet, >1 threshold = Doublet/Multiplet.
5. Remove unclassified cells, repeat steps 1-4 until all cells are classified as either singlets or doublets/multiplets.

![alternativetext](/Figures/MULTIseq_ClassificationWorkflow.png)

Notably, assigning local maxima to background and positive peaks assumes that the background peak will always be the mode. As one can imagine, this assumption breaks down when (i) you have only 2 samples and/or (ii) one sample is more frequent than the sum of all other sample frequencies. If you believe this is happening in your dataset, I advise either (i) using one of the alternative sample classification algorithms linked below or (ii) manually thresholding your barcode distributions, tracking cells surpassing each manually-set thresheold, and classifying singlets/doublets/negatives as described in step 4.

## MULTIseq semi-supervised negative-cell reclassification 
In its current form, MULTI-seq barcoding is an imperfect process that produces a small fraction of cells that cannot be classified into sample groups. These negative cells are of two varieties: True and false negatives. True negatives result from cells with poor barcode labeling. In contrast, false negatives result from algorithmic misclassification. Since a single inter-maxima quantile threshold is applied to all barcodes during sample classification, we believe false negatives arise because this thresholding strategy may be sub-optimal for a subset of barcode distributions. Thus, although false negatives have poor *absolute* signal in comparison to high-confidence singlets, we reasoned that false negatives could be ‘rescued’ by computing the *relative* strength of each barcode signal on a cell-by-cell basis.

Broadly, negative cell reclassification uses the initial MULTI-seq sample classification results as 'ground-truth' labels during semi-supervised k-means clustering of negative cells (as in Cell Hashing). This pipeline can be split into seven distinct steps:  
1. Record the total number of thresholds that each negative cell surpasses at each inter-maxima quantile.
2. Compute each negative cell’s classification stability (CS) – i.e., the number of quantiles across which a cell surpasses a single threshold.
3. Subset equal numbers of ‘ground-truth’ cells from the original classification results.
4. Perform semi-supervised k-means clustering on merged data including ‘ground-truth’ and negative cells. Clustering is semi-supervised because one member of each ‘ground-truth’ sample group is used to initialize cluster centers.
5. Compute the rate at which ‘ground-truth’ and negative cell classifications match the k-means results.
6. Iteratively repeat steps 4 and 5 using a different ‘ground-truth’ cell to initialize cluster centers during each iteration. Repeat until all ‘ground-truth’ cells have been used.
7. Compare k-means matching rates between ‘ground-truth’ and negative cells binned according to CS values. Negative cells with CS values resulting in matching rates that approximate ‘ground-truth’ matching rates are reclassified.

![alternativetext](/Figures/MULTIseq_NegativeCellReclassification.png)

# Tutorial: 96-plex HMEC sample multiplexed scRNA-seq
## Step 1: MULTI-seq sample barcode pre-processing and alignment
```R
## Define vectors for reference barcode sequences and cell IDs
bar.ref <- load("/path/to/BClist.Robj")
cell.id.vec <- load("/path/to/cell.id.vec.Robj")
```
![alternativetext](/Figures/Tutorial_BarRef_CellIDs.png)

```R
## Pre-process MULTI-seq sample barcode FASTQs
readTable <- MULTIseq.preProcess(R1 = '/path/to/R1.fastq.gz', R2 = '/path/to/R2.fastq.gz', cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8))
```
![alternativetext](/Figures/Tutorial_preprocess.out.png)

```R
## Perform MULTI-seq sample barcode alignment
bar.table <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
```
![alternativetext](/Figures/Tutorial_align.out.png)

## Step 2: Visually inspect sample barcode quality
Prior to sample classification, it is important to check that your data includes cells that were successfully labeled with every intended barcode. Missing barcodes will severely impair sample classification performance. Thus, we suggest visually inspecting whether every barcode is enriched in a unique region of barcode space (i.e., as visualized with tSNE).

```R
## Visualize barcode space
bar.tsne <- barTSNE(bar.table[,1:96]) 
## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell. 

pdf("bc.check.pdf")
for (i in 3:ncol(bar.tsne)) {
    g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
    print(g)
}
dev.off()
```
![alternativetext](/Figures/Tutorial_good.vs.bad.bc.tsnes.png)

Remove columns in 'bar.table' corresponding to missing barcodes prior to beginning sample classification. 

## Step 3: Sample Classification
```R
## Round 1 -----------------------------------------------------------------------------------------------------
## Perform Quantile Sweep
bar.table.full <- bar.table[,1:96]
good.bars <- paste("Bar",1:90,sep="")  # NOTE: In this hypothetical example, barcodes 91-96 were not detected
bar.table <- bar.table.full[, good.bars]  # Remove missing bars and summary columns
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
```
![alternativetext](/Figures/Tutorial_good.qsweep.png)

If this plot does not resemble the visualized distribution above, there may be missing barcodes remaining.

```R
## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round2.calls),neg.cells)

```
![alternativetext](/Figures/Tutorial_final.classification.results.png)

## Step 4 (optional): Semi-Supervised Negative Cell Reclassification
```R
## Perform semi-supervised negative cell reclassification
reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)
```
![alternativetext](/Figures/Tutorial_rescue.out.png)

```R
## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
    geom_point() + xlim(c(nrow(reclass.res)-1,1)) + 
    ylim(c(0,1.05)) +
    geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)
```
![alternativetext](/Figures/Tutorial_match.rate.cs.png)

```R
## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 16) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
```

## Other Sample Demultiplexing Methods
[CITE-Seq Count](https://github.com/Hoohm/CITE-seq-Count)
[Seurat](https://satijalab.org/seurat/v3.1/hashing_vignette.html)
[DemuxEM](https://pegasus.readthedocs.io/en/0.x/usage.html?highlight=demuxEM#pegasus-demuxem)
[Solo](https://github.com/calico/solo)
[GMM-Demux](https://github.com/CHPGenetics/GMM-demux)

## Referencens
1.  Stoeckius M, Zheng S, Houck-Loomis B, Hao S, Yeung BZ, Smibert P, Satija R. Cell "hashing" with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. Genome Biology. 2018. 19(1):224.

2.  Adamson B, Norman TM, Jost M, Cho MY, Nuñez JK, Chen Y, et al. A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response. Cell. 2016. 167(7):1867-82.e21.

3.  Dixit A, Parnas O, Li B, Chen J, Fulco CP, Jerby-Arnon L, et al. Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens. Cell. 2016. 167(7):1853-66.e17.

4.  Gaublomme JT, Li B, McCabe C, Knecht A, Yang Y, Drokhlyansky E, Van Wittenberghe N, Waldman J, Dionne D, Nguyen L, De Jager PL, Yeung B, Zhao X, Habib N, Rozenblatt-Rosen O, Regev A. Nuclei multiplexing with barcoded antibodies for single-nucleus genomics. 2019. Nature Communications. 10(1):2907.

5.  Bernstein NJ, Fong NL, Lam I, Roy MA, Hendrickson DG, Kelley DR. Solo: Doublet Identification in Single-Cell RNA-Seq via Semi-Supervised Deep Learning. Cell Systems. 2020. S2405-4712(20)30195-2.

6.  Gehring J, Park JH, Chen S, Thomson M, Pachter L. Highly multiplexed single-cell RNA-seq by DNA oligonucleotide tagging of cellular proteins. Nature Biotechnology. 2019. 38(1):35-8.
