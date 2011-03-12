\name{mRNAdoubleKDefficiency}
\alias{mRNAdoubleKDefficiency}
\docType{data}
\title{
mRNA levels for double knock downs
}
\description{
qPCR measurements for the mRNA level after a double gene knock down (ratio relative to wild type control). The experiments tests the knock down efficiency in the presence of a second gene knock down.
}
\usage{data(mRNAdoubleKDefficiency)}
\format{
  A data frame with 320 observations on the following 5 variables.
  \describe{
    \item{\code{template}}{an ordered factor with levels \code{Fluc} < \code{CG10417} < \code{CG13197} < \code{CG9391} < \code{egr} < \code{lic} < \code{PRL-1} < \code{Rho1} < \code{Tak1}}
    \item{\code{query}}{a factor with levels \code{CG10417 [query dsRNA]} \code{CG13197 [query dsRNA]} \code{CG9391 [query dsRNA]} \code{egr [query dsRNA]} \code{lic [query dsRNA]} \code{PRL-1 [query dsRNA]} \code{Rho1 [query dsRNA]} \code{Tak1 [query dsRNA]}}
    \item{\code{qPCR.target}}{a factor with levels \code{CG10417} \code{CG13197} \code{CG9391} \code{egr} \code{lic} \code{PRL-1} \code{Rho1} \code{Tak1}}
    \item{\code{passage}}{a factor with levels \code{passage 4} \code{passage 42}}
    \item{\code{RNAi}}{a numeric vector}
  }
}
\source{
Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, Figure S2.
}
\examples{
data(mRNAdoubleKDefficiency)
head(mRNAdoubleKDefficiency)
}
\keyword{datasets}
