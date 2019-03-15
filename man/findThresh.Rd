\name{ findThresh }

\alias{ findThresh }

\title{ findThresh }

\description{ 'findThresh' summarizes the results of the quantile sweep performed during each round of the MULTI-seq sample classification workflow. Output is used to select the barcode-specific threshold used to classify cells during each round. }

\usage{ findThresh(call.list) }

\arguments{
  \item{call.list}{ List containing the singlet/doublet/negative classifications for each quantile threshold tested during a round of the MULTI-seq sample classification workflow.  }
}

\details{ Requires the 'reshape2' R package. }

\value{ Returns a list containing (1) 'res': proportion of cells called as singlets (pSinglet), doublets (pDoublet), and negatives (pNegative) at each quantile tested, and (2) 'extrema': the local maxima for the pSinglet qunatile distribution. }

\references{ Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. }

\author{ Chris McGinnis }

\examples{
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by=0.02)) {
    print(q)
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(data=barTable, q=q)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
  }

  threshold.results1 <- findThresh(bar.table_sweep.list)
  round1.calls <- classifyCells(barTable, q=findQ(threshold.results1$res, threshold.results1$extrema))
}
