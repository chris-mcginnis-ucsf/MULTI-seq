\name{ localMaxima }

\alias{ localMaxima }

\title{ localMaxima }

\description{ 'localMaxima' is an internal function to 'classifyCells' and 'findReclassCells' which identifies local maxima amongst sample barcode UMI probability density functions. }

\usage{ localMaxima(x) }

\arguments{
  \item{x}{ Probability-density function, as generated using 'bkde'. }
}

\details{ As implemented by 'Tommy' at https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima }

\value{ Vector of local maxima }

\author{ Chris McGinnis }

\examples{
  model <- approxfun(bkde(barTable.n[,i], kernel="normal"))
  x <-  seq(from=quantile(barTable.n[,i],0.001), to=quantile(barTable.n[,i],0.999),
    length.out=100)
  extrema <- localMaxima(model(x))
}

