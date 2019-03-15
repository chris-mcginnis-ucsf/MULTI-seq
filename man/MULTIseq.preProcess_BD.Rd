\name{ MULTIseq.preProcess_BD }

\alias{ MULTIseq.preProcess_BD }

\title{ MULTIseq.preProcess_BD }

\description{ 'MULTIseq.preProcess_BD' reads in raw BD HashTag sample barcode FASTQs and allocates reads into cell barcode, UMI, and HashTag subsets. \cr

R1 reads must exactly match cell IDs provided in whitelist. Specifies cell barcode from R1 reads 1-16, UMI from R1 reads 17-26/28 (V2/V3 chemistry), and sample barcode from R2 reads 1-70 (includes 25bp constant and 45bp HashTag sequences.
}

\usage{ MULTIseq.preProcess_BD(R1, R2, cellIDs, chemistry="V3") }

\arguments{
  \item{R1} { Read 1 raw FASTQ file path }
  \item{R2} { Read 2 raw FASTQ file path }
  \item{cellIDs} { Vector of cell IDs used to parse R1 FASTQ file }
  \item{chemistry} { 10X reagent chemistry, determines UMI sequence position. Set to either "V2" or "V3" }
}

\details{ Requires 'ShortRead' R package. }

\value{ An nRead x 3 dataframe with columns correpsonding to (1) cell barcode, (2) UMI, and (3) HashTag constant and sample barcode sequence. }

\references{ Morgan M, Anders S, Lawrence M, Aboyoun P, Pagès H, Gentleman R. ShortRead: a Bioconductor package for input, quality assessment and exploration of high-throughput sequence data. Bioinformatics. 2009; 25:2607-8. }

\author{ Chris McGinnis}

\examples{ readTable <- MULTIseq.preProcess_BD(R1 = '/path/to/R1.fastq.gz', R2 = '/path/to/R2.fastq.gz', cellIDs = cell.id.vec, chemistry = "V3") }
