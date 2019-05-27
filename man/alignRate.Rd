\name{alignRate}

\alias{alignRate}

\title{alignRate}

\description{ 'alignRate' computes the proportion of sample barcode reads that align with <1 mismatch to any sample barcode reference sequence. }

\usage{ alignRate(readTable, cellIDs, ref) }

\arguments{
  \item{readTable}{ Read table comprised of CellID (R1), UMI (R1), and sample barcode (R2), as generated using MULTIseq.preProcess.}
  \item{cellIDs}{ Vector of cellIDs with which to compute the alignment rate. }
  \item{ref}{ Vector of sample barocde sequences with which to compute the alignment rate. }
}

\value{ Vector of length cellIDs with alignment rate values (0-100) for each cell. }

\author{ Chris McGinnis }

\examples{ ar <- alignRate(readTable, cellIDs, ref) }

