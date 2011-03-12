\name{mRNAsingleKDefficiency}
\alias{mRNAsingleKDefficiency}
\docType{data}
\title{
mRNA levels for single gene knock downs
}
\description{
qPCR measurements for the mRNA level after a single gene knock down (ratio relative to wild type control). The experiments is done for two independent designs of RNAi reagents.
}
\usage{data(mRNAsingleKDefficiency)}
\format{
  A data frame with 89 observations on the following 5 variables.
  \describe{
    \item{\code{Symbol}}{a character vector}
    \item{\code{MeanDesign1}}{a numeric vector}
    \item{\code{StderrDesign1}}{a numeric vector}
    \item{\code{MeanDesign2}}{a numeric vector}
    \item{\code{StderrDesign2}}{a numeric vector}
  }
}
\source{
Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, Figure S1.
}
\examples{
data(mRNAsingleKDefficiency)
head(mRNAsingleKDefficiency)
}
\keyword{datasets}
