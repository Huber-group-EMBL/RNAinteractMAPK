\name{Networks}
\alias{Networks}
\docType{data}
\title{
 Knock interaction networks
}
\description{
This dataset is a subset of the DroID database. It contains the known (genetic) interactions between the genes regarded in the main screen.
}
\usage{data(Networks)}
\format{
  A data frame with 402 observations on the following 5 variables.
  \describe{
    \item{\code{gene1}}{a character vector}
    \item{\code{gene2}}{a character vector}
    \item{\code{correlation}}{a numeric vector}
    \item{\code{genetic}}{a numeric vector}
    \item{\code{human}}{a numeric vector}
  }
}
\source{
Data as used in Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, Figure S12.

The data is a subset from the drosophila interactions database (DroID), http://www.droidb.org, Data version 2010_10 updated 20 October 2010.

Murali T, Pacifico S, Yu J, Guest S, Roberts GG 3rd, Finley RL Jr . DroID 2011: a comprehensive, integrated resource for protein, transcription factor, RNA and gene interactions for Drosophila. Nucleic Acids Res. 2010 Oct 29.
}
\examples{
data(Networks)
head(Networks)
}
\keyword{datasets}
