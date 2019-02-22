\name{ findThresh }

\alias{ findThresh }

\title{ findThresh }

\description{ 'findThresh' summarizes the results of the quantile sweep performed during each round of the MULTI-seq sample classification workflow. Output is used to select the barcode-specific threshold used to classify cells during each round. }

\usage{ findThresh(call.list, id) }

\arguments{
  \item{call.list}{ List of Singlet/Doublet/Negative classifications for each quantile threshold tested during a round of the MULTI-seq sample classification workflow. }
  \item{id}{ String appended to 'findThresh' outputs. Usually corresponds to sample classification round (e.g., "round1", "round2", "round3", etc.) }
}

\details{ Requires the 'reshape2' R package. }

\value{ Sets two objects to the global environment: res_id and extrema_id, where 'id' is set to the 'id' argument (e.g., id="round1 results in res_round1 and extrema_round1 assignments). 'res_id' is a dataframe describing pSinglet, pNegative, and pDoublet across all barcodes at every quantile tested. 'extrema_id' are the local pSinglet maxima which, when maximized, are used to define barcode-specific thresholds prior to sample classification. }

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

  findThresh(bar.table_sweep.list, "round1") ## sets res_round1, extrema_round1 to global env
  round1.calls <- classifyCells(barTable, q=findQ(res_round1, extrema_round1))
}
