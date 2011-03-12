\name{MAPK.estimate.TPS}
\alias{MAPK.estimate.TPS}
\alias{MAPK.cv.TPS}
\alias{MAPK.plot.TPS.all}
\alias{MAPK.plot.TPS.single}
\alias{MAPK.screen.as.array}
\title{
genetic interaction surfaces
}
\description{
Genetic interaction surfaces are estimated from a dilution experiment. Cells are treated with two RNAi's. The concentration of the RNAi reagent is changed in 8 steps. All 8 x 8 combinations of concentrations are tested for 6 x 6 gene pairs.  
}
\usage{
MAPK.screen.as.array(data, Anno)
MAPK.estimate.TPS   (A, DF, n.out = 8, channel = 1)
MAPK.cv.TPS         (A, range.df = 6:56, channel = 1)
MAPK.plot.TPS.all   (TPSmodel, range = c(-6, 6), fill = c("cornflowerblue",
                     "cornflowerblue", "black", "#777700", "yellow"),
		     channel = 1)
MAPK.plot.TPS.single(gene1, gene2, TPSmodel, range = c(-6, 6),
                     fill = c("cornflowerblue", "cornflowerblue", "black",
		     "#777700", "yellow"), channel = 1)
}
\arguments{
data, Anno
  \item{data}{A data.frame containing the read.out of the dilution screen. Each row is one well. Each column one feature.}
  \item{Anno}{A data.frame containing the plate configuration. For each row in data there should a row in Anno.}
  \item{A}{An array of dimension concentration x concentration x genes x genes x channel as returned by \code{MAPK.screen.as.array}.}
  \item{DF}{A 6 x 6 matrix of degrees of freedom for the thin spline plate regression.}
  \item{n.out}{number of points for sampling from the regression function.}
  \item{channel}{The channel used.}
  \item{range.df}{The range of degrees of freedom that is considered for cross validation.}
  \item{range}{The range of pairwise interaction scores that is shown by the colorbar.}
  \item{gene1, gene2}{The genes for which the interaction surface is plotted.}
  \item{TPSmodel}{The TPS model estimated by \code{MAPK.estimate.TPS}}
  \item{fill}{The range of colors used for the color code.}
}
\details{
The screen readout can be reshaped as an array with dimensions concentration x concentration x genes x genes x channel by the function \code{MAPK.screen.as.array}. Then the function \code{MAPK.estimate.TPS} fits a regression in the 8 x 8 pairwise dilution series. The degrees of freedom for the regression can be estimated automatically by cross validation with the function \code{MAPK.cv.TPS}. Finally one can plot the interaction surface for a single gene or an overview of interaction surfaces for all genes with the functions \code{MAPK.plot.TPS.single} or \code{MAPK.plot.TPS.all}.
}
\value{
\itemize{
\item \code{MAPK.screen.as.array} returns an array of dimensions concentration x concentration x genes x genes x channel with the screen data.

\item \code{MAPK.estimate.TPS} returns a regression model estimated by thin plate splines for each pair of genes and subsampled matrices.

\item \code{MAPK.cv.TPS.all} returns 
\itemize{
\item \code{DF} a matrix with degrees of freedom.
\item \code{CVerror} The prediction error estimated by cross validation.
\item \code{CVerrorSD} The standard deviation of the prediction error estimated by cross validation.
}

\item \code{MAPK.plot.TPS.single},\code{MAPK.plot.TPS.all}: An object of class "trellis". See \code{\link{levelplot}} for details.
}
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}, \code{\link{RNAinteractMAPK-package}}
}
\keyword{ misc }
