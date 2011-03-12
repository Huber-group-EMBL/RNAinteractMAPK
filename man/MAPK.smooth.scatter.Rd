\name{MAPK.smooth.scatter}
\alias{MAPK.smooth.scatter}
\title{
smooth scatter using grid raster
}
\description{
This function is a reimplementation of \code{\link{smoothScatter}}. For nicer graphics output the background image is written by \code{\link{grid.raster}}. It is recommended to use the \code{\link{smoothScatter}} function from the \code{\link{graphics}} package.
}
\usage{
MAPK.smooth.scatter(x, y, n = 75,
                    nrpoints = 100, col = "blue", 
		    pch = 20, size = unit(0.3, "char"), cex = 1.2,
		    colramp = colorRampPalette(c("white", "blue", "green", "yellow", "red"))(256),
		    xlab = "", ylab = "", respect = FALSE)
}
\arguments{
  \item{x}{x-values.}
  \item{y}{y-values. Has to be the same length as x.}
  \item{n}{nr of bins used for the kernel density estimation.}
  \item{nrpoints}{nrpoints points in the lowest density region will be plotted. This allows the identification of outliers.}
  \item{col}{color of points.}
  \item{pch}{symbol to plot points.}
  \item{size}{The size of the points.}
  \item{cex}{The size of the label text.}
  \item{colramp}{color ramp for the density plot.}
  \item{xlab,ylab}{axis labels.}
  \item{respect}{A logical value indicating if the height and width of the axis scales should respect each other.}
}
\details{
Plots a density plot with \code{\link{grid}} graphics.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteractMAPK-package}}, \code{\link{smoothScatter}}
}
\keyword{ misc }
