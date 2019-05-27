\name{ MULTIseq.preProcess }

\alias{ MULTIseq.preProcess }

\title{ MULTIseq.preProcess }

\description{ 'MULTIseq.preProcess' reads in raw MULTI-seq sample barcode FASTQs and allocates reads into cell barcode, UMI, and sample barcode subsets. }

\usage{ MULTIseq.preProcess(R1, R2, cellIDs, cell=c(1,16), umi=c(17,28), tag=c(1,8)) }

\arguments{
  \item{R1}{ Read 1 raw FASTQ file path }
  \item{R2}{ Read 2 raw FASTQ file path }
  \item{cellIDs}{ Vector of cell IDs used to parse R1 FASTQ file }
  \item{cell}{ Numerical vector of length 2 specifying beginning and end position of cell barcode in R1 (Default = 1:16)}
  \item{umi}{ Numerical vector of length 2 specifying beginning and end position of UMI in R1 (Default = 17:28) }
  \item{tag}{ Numerical vector of length 2 specifying beginning and end position of sample tag in R2 (Default = 1:8) }
}

\details{ Requires 'ShortRead' R package. }

\value{ An nRead x 3 dataframe with columns correpsonding to (1) cell barcode, (2) UMI, and (3) MULTI-seq sample barcode sequences. }

\references{ Morgan M, Anders S, Lawrence M, Aboyoun P, Pagès H, Gentleman R. ShortRead: a Bioconductor package for input, quality assessment and exploration of high-throughput sequence data. Bioinformatics. 2009; 25:2607-8. }

\author{ Chris McGinnis }

\examples{ readTable <- MULTIseq.preProcess(R1 = '/path/to/R1.fastq.gz', R2 = '/path/to/R2.fastq.gz', cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8)) }

