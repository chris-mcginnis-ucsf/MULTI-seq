\name{ MULTIseq.preProcess }

\alias{ MULTIseq.preProcess }

\title{ MULTIseq.preProcess }

\description{ 'MULTIseq.preProcess' reads in raw MULTI-seq sample barcode FASTQs and allocates reads into cell barcode, UMI, and sample barcode subsets. \cr

R1 reads must exactly match cell IDs provided in whitelist. Specifies cell barcode from R1 reads 1-16, UMI from R1 reads 17-26/28 (V2/V3 chemistry), and sample barcode from R2 reads 1-8.
}

\usage{ MULTIseq.preProcess(R1, R2, cellIDs, chemistry="V3") }

\arguments{
  \item{R1} { Read 1 raw FASTQ file path }
  \item{R2} { Read 2 raw FASTQ file path }
  \item{cellIDs} { Vector of cell IDs used to parse R1 FASTQ file }
  \item{chemistry} { 10X reagent chemistry, determines UMI sequence position. Set to either "V2" or "V3" }
}

\details{ Requires 'ShortRead' R package. }

\value{ An nRead x 3 dataframe with columns correpsonding to (1) cell barcode, (2) UMI, and (3) MULTI-seq sample barcode sequences. }

\references{ Morgan M, Anders S, Lawrence M, Aboyoun P, Pagès H, Gentleman R. ShortRead: a Bioconductor package for input, quality assessment and exploration of high-throughput sequence data. Bioinformatics. 2009; 25:2607-8. }

\author{ Chris McGinnis }

\examples{ readTable <- MULTIseq.preProcess(R1 = '/path/to/R1.fastq.gz', R2 = '/path/to/R2.fastq.gz', cellIDs = cell.id.vec, chemistry = "V3") }

