\name{PhysicalInteractions}
\alias{PhysicalInteractions}
\docType{data}
\title{
Known physical interactions
}
\description{
This dataset contains a collection of known physical interactions assembled from the literature. It contains the known pathway structure of the RasMAPK and the JNK pathway.
}
\usage{data(PhysicalInteractions)}
\format{
  A data frame with 29 observations on the following 2 variables.
  \describe{
    \item{\code{V1}}{a character vector}
    \item{\code{V2}}{a character vector}
  }
}
\source{
Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, Figure S13.
}
\examples{
data(PhysicalInteractions)
head(PhysicalInteractions)
}
\keyword{datasets}
