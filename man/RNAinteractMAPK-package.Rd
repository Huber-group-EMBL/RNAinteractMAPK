\name{RNAinteractMAPK-package}

\alias{RNAinteractMAPK-package}
\alias{RNAinteractMAPK}

\docType{package}

\title{Mapping of Signalling Networks through Synthetic Genetic Interaction Analysis by RNAi.}

\description{
The package contains the data and the source code to reproduce the results and figures from the paper

T. Horn, T. Sandmann, B. Fischer, W. Huber, M. Boutros. Mapping of Signalling Networks through Synthetic 
Genetic Interaction Analysis by RNAi. Nature Methods, 2011.
}

\section{Package content}{

  See vignette("RNAinteractMAPK") for more detail on how to obtain the data used for specific figures. In addition this vignette contains the complete analysis and the generation of all figures. 

  The main screen can be loaded by \code{data("Dmel2PPMAPK", package="RNAinteractMAPK")}. Access to the pairwise interaction data is done via the \code{\link{getData}} function from the \code{\link{RNAinteract}}-package. See example below.

  The following datasets are provided with this package:
  \tabular{ll}{
	\code{\link{Dmel2PPMAPK}} \tab interaction data of main screen. See example below. \cr
	\code{\link{ElpB1phenotype}} \tab in vivo experiment on ectopic wing vein formation (used in Figure 5f)\cr
	\code{\link{mRNAdoubleKDefficiency}} \tab mRNA level after double gene knockdown (used in Figure S1) \cr
	\code{\link{mRNAsingleKDefficiency}} \tab mRNA level after single gene knockdown (used in Figure S2) \cr
	\code{\link{singleKDphenotype}} \tab single knockdown phenotypes (used in Figure S3) \cr
	\code{\link{dsRNAiDilutionSeries}} \tab dsRNA dilution series (used in Figure S5) \cr
  	\code{\link{Networks}} \tab known interactions from DroID (used in Figure S12)  \cr
	\code{\link{pathwayMembership}} \tab pathway membership of genes (used in Figure S13) \cr
	\code{\link{PhysicalInteractions}} \tab Known physical interactions (used in Figure S13)  \cr
	\code{\link{cellTiterGlo}} \tab cellTiterGlo viability data (used in Figure S15)
  }


  Within this package a number of specialized functions is defined written for the analysis of the MAPK interaction screen and additional experiments shown in the paper. These functions are not intended to be general purpose analysis functions, but should document the steps of analysis of the paper. Therefore, these functions are not exported. A list of functions is given below.
  A general purpose package for the analysis of genetic interaction screens is the package \code{\link{RNAinteract}}.

  The following functions are provided within this package. 

  Functions used for the classification: \code{MAPK.predict.classification}, \code{MAPK.cv.classifier}, \code{MAPK.getCV}, \code{MAPK.ternary.plot}, \code{MAPK.getXY}, \code{MAPK.plot.classification}.

  Functions for the analysis of the interaction surfaces: \code{MAPK.plot.TPS.single}, \code{MAPK.plot.TPS.all}, \code{MAPK.estimate.TPS}, \code{MAPK.cv.TPS}, \code{MAPK.screen.as.array}.
  
  A function to write the hit list: \code{MAPK.report.gene.lists.paper}.

  A function to plot a heatmap: \code{MAPK.plot.heatmap.raster}.

  A function to plot smaooth scatters: \code{MAPK.smooth.scatter}.

}

\details{
See vignette("RNAinteractMAPK") for details.
}

\author{
Bernd Fischer

Maintainer: Bernd Fischer <bernd.fischer@embl.de>
}
\references{
T. Horn, T. Sandmann, B. Fischer, W. Huber, M. Boutros. Mapping of Signalling Networks through Synthetic 
Genetic Interaction Analysis by RNAi. Nature Methods, 2011.
}
\keyword{ package
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
~~ documentation directory ~~
 }
\seealso{
\code{\link{RNAinteractMAPK-package}}
}
\examples{
data(Dmel2PPMAPK)
Dmel2PPMAPK

# Obtain the pairwise interaction matrix
PI <- getData(Dmel2PPMAPK, type="pi", format="targetMatrix", screen="mean", withoutgroups = c("pos", "neg"))
}



