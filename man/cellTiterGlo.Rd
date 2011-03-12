\name{cellTiterGlo}
\alias{cellTiterGlo}
\docType{data}
\title{
Comparison of interaction experiment with an cellTiterGlo viability assay
}
\description{
For ten gene pairs genetic interactions are measured. The experiment contains 24 different conditions. These are repeated in each row of the three 384 mutli well plates. The data.frame contains the plate annotation as well as the viability readout for the three plates.
}
\usage{data(cellTiterGlo)}
\format{
  A data frame with 384 observations on the following 6 variables.
  \describe{
    \item{\code{Well}}{a character vector}
    \item{\code{dsRNA_1}}{a character vector}
    \item{\code{dsRNA_2}}{a character vector}
    \item{\code{plate1}}{a numeric vector}
    \item{\code{plate2}}{a numeric vector}
    \item{\code{plate3}}{a numeric vector}
  }
}
\source{
Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, Figure S15.
}
\examples{
data(cellTiterGlo)
head(cellTiterGlo)
}
\keyword{datasets}
