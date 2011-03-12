\name{singleKDphenotype}
\alias{singleKDphenotype}
\alias{singleKDphenotypeAnno}
\docType{data}
\title{
Single knockdown phenotype
}
\description{
This data.frame \code{singleKDphenotype} contains a screen assessing the single knock down phenotypes (nrCells, intensity, and area) of the tested genes. \code{singleKDphenotypeAnno} is a data.frame describing the plate annotation.
}
\usage{data(singleKDphenotype)}
\format{
  The format is:
 chr "singleKDphenotype"
 chr "singleKDphenotypeAnno"
}
\source{
Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, Figure S3.
}
\examples{
data(singleKDphenotype)
head(singleKDphenotype)
head(singleKDphenotypeAnno)
}
\keyword{datasets}
