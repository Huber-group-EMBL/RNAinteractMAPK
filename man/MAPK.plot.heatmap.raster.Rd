\name{MAPK.plot.heatmap.raster}
\alias{MAPK.plot.heatmap.raster}
\title{
Plots a heatmap using grid.raster
}
\description{
This functions provides a grid plot that displays the raster image of a heatmap without any axis or label. This function is adapted from the function \code{\link{grid.sgiHeatmap}} from the package \code{\link{RNAinteract-package}}. It is highly recommended to use the original function \code{\link{grid.sgiHeatmap}}.
}
\usage{
MAPK.plot.heatmap.raster(X, subset = NULL,
                         hc.row = NULL, hc.col = NULL,
                         pi.max = NULL)
}
\arguments{
  \item{X}{A matrix of pairwise interaction scores.}
  \item{subset}{A subset of genes that are displayed in the rows.}
  \item{hc.row}{A hclust object.}
  \item{hc.col}{A hclust object.}
  \item{pi.max}{The maximum interaction score of the colorbar. All interaction scores larger than this value will be displayed in the same color.}
}
\value{Nothing is returned.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}, \code{\link{RNAinteractMAPK-package}}
}
\keyword{ misc }
