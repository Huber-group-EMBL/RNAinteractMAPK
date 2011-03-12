\name{dsRNAiDilutionSeries}
\alias{dsRNAiDilutionSeries}
\alias{dsRNAiDilutionSeriesAnno}
\alias{dsRNAiDilutionSeriesDF}
\docType{data}
\title{
dsRNA dilution series
}
\description{
A dilution series for 6 x 6 gene. For each gene pair all combinations of 8 different concentrations of dsRNA reagent are screened. Three readout channels (nrCells, area, intensity) are available in the \code{data.frame dsRNAiDilutionSeries}. The plate annotation is given in the \code{data.frame dsRNAiDilutionSeriesAnno} and precomputed degrees of freedom for thin plate splines are available in the \code{matrix dsRNAiDilutionSeriesDF}.
}
\usage{data(dsRNAiDilutionSeries)}
\format{
  \itemize{
  \item \code{dsRNAiDilutionSeries}:
    A data frame with 2688 observations on the following 3 variables.
    \describe{
      \item{\code{nrCells}}{a numeric vector representing the number of cells readout.}
      \item{\code{area}}{a numeric vector representing the are readout.}
      \item{\code{intensity}}{a numeric vector representing the intensity readout.}
    }
  \item \code{dsRNAiDilutionSeriesAnno}:
    A data frame with 2688 observations on the following 7 variables.
    \describe{
      \item{\code{plate}}{a numeric vector representing the plate number.}
      \item{\code{well}}{a numeric vector representing the well.}
      \item{\code{row}}{a numeric vector representing the row on plate.}
      \item{\code{col}}{a numeric vector representing the column on plate.}
      \item{\code{targetID1}}{a numeric vector representing the target number of the first reagent (see targets).}
      \item{\code{targetID2}}{a numeric vector representing the target number of the second reagent (see targets).}
      \item{\code{targets}}{a data.frame representing the 49 target reagents For each target reagent the gene name and the dsRNAi concentration is given.}
    }
  \item \code{dsRNAiDilutionSeriesDF}:
    A 6 x 6 matrix with degrees of freedom for thin-plate spline regression.
  }
}
\source{
Horn, Sandmann, Fischer, Huber, Boutros, Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi, Nature Methods, 2011, Figure 1abc and Figure S5.
}
\examples{
data(dsRNAiDilutionSeries)
head(dsRNAiDilutionSeries)
head(dsRNAiDilutionSeriesAnno)
head(dsRNAiDilutionSeriesDF)
}
\keyword{datasets}
