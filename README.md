# MULTI-seq Sample Classification
MULTI-seq is a method for single-cell RNA sequencing sample multiplexing (for more information, check out our preprint: https://www.biorxiv.org/content/early/2018/08/08/387241). 

Here, we provide source code for our sample classification framework. Our workflow expands on concepts borrowed from Cell-Hashing (Stoeckius et al., 2017) and Perturb-seq (Adamson et al., 2016; Dixit et al., 2016). As in Perturb-Seq, we model bimodal barcode UMI distributions across all cells using Guassian-kernel density estimation. We then define and optimize barcode-specific thresholds according to these models. Then, as in Cell-Hashing, we classify cells according to which threshold it surpassess, and define doublets as cells surpassing >1 treshold.
  
# How to run
1. Generate barcode count matrix (e.g., using CITE-seq Count: https://github.com/Hoohm/CITE-seq-Count) 
2. Download and unzip 'MULTIseq.zip'
3. Load 'SampleClassificationScript.R' into R -- Follow the instructions
 
![alternativetext](Workflow.png)

# Referencens
1. Stoeckius M, Zheng S, Houck-Loomis B, Hao S, Yeung BZ, Smibert P, Satija R. Cell "hashing" with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. 2017. Preprint. bioRxiv doi: 10.1101/237693.
2. Adamson B, Norman TM, Jost M, Cho MY, Nuñez JK, Chen Y, et al. A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response. Cell. 2016; 167(7):1867-82.e21.
3. Dixit A, Parnas O, Li B, Chen J, Fulco CP, Jerby-Arnon L, et al. Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens. Cell. 2016; 167(7):1853-66.e17.
