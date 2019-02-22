\name{ bucket_readTable }

\alias{ bucket_readTable }

\title{ bucket_readTable }

\description{ 'bucket_readTable' is an internal function for 'MULTIseq.align' which splits read tables into 10,000-cellID buckets. Requires input from 'bucket_cellIDs' }

\usage{ bucket_readTable(readTable, cellIDs_bucketed) }

\arguments{
  \item{readTable}{ Parsed read table, as provided by 'MULTIseq.preProcess' }
  \item{cellIDs_bucketed}{ Bucketed list of cellIDs, as provided by 'bucket_cellIDs' }
}

\value{ A list of readTables containing all reads associated with sets of 10,000-cells. }

\author{ Chris McGinnis }

\examples{ readTable.list <- bucket_readTable(readTable, cellIDs_bucketed) }
