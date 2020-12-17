\name{Dmel2PPMAPK}
\alias{Dmel2PPMAPK}
\docType{data}
\title{
The interaction data of the main screen
}
\description{
\code{Dmel2PPMAPK} is an object of class \code{\link{RNAinteract}}. It contains the raw data, the computed main effects, pairwise interaction scores, p-values and q-values computed by a t-test. The package vignette vontains the complete code and documentation for the statistical analysis.
}
\usage{data(Dmel2PPMAPK)}
\format{
An object of class \code{\link{RNAinteract}}.
}
\source{
Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, main interaction screen.
}
\examples{
data(Dmel2PPMAPK)
Dmel2PPMAPK

# Obtain the pairwise interaction matrix
PI <- getData(Dmel2PPMAPK, type="pi", format="targetMatrix", 
              screen="mean", withoutgroups = c("pos", "neg"))
}
\keyword{datasets}
