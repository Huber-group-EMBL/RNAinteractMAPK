\name{MAPK.report.gene.lists.paper}
\alias{MAPK.report.gene.lists.paper}
\title{
A hitlist report.
}
\description{
Reports the hitlist of genetic interactions, with p-values from a t-test with pooled variance estimate, from limma, and from Hotelling T^2 test.
}
\usage{
MAPK.report.gene.lists.paper(sgi, sgilimma, sgi3T2, screen = "mean")
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}} containing p-values from a t-test with pooled variance estimate.}
  \item{sgilimma}{An object of class \code{\link{RNAinteract}} containing p-values from limma.}
  \item{sgi3T2}{An object of class \code{\link{RNAinteract}} containing p-values from a Hotelling T^2 test.}
  \item{screen}{The screen name for which the report should be written.}
}
\details{
Writes tab-separated lists for each single test as well as a joint table with all three tests.
}
\value{
Nothing is returned.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}, \code{\link{RNAinteractMAPK-package}}
}
\keyword{ misc }
