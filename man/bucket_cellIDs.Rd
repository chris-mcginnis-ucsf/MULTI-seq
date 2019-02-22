\name{bucket_cellIDs}

\alias{bucket_cellIDs}

\title{bucket_cellIDs }

\description{ 'bucket_cellIDs' is an internal function for 'MULTIseq.align' which splits cellIDs into 10,000-cell buckets which are subsequently used to parse readTables. }

\usage{ bucket_cellIDs(cellIDs) }

\arguments{
  \item{cellIDs}{ Vector of cell IDs to be bucketed. }
}

\value{ A list of cellID vectors of length 10,000. }

\author{ Chris McGinnis }

\examples{ cellID_list <- bucket_cellIDs(cellIDs) }
