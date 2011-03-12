\name{MAPK.cv.classifier}
\alias{MAPK.cv.classifier}
\alias{MAPK.predict.classification}
\alias{MAPK.plot.classification}
\title{
A classifier for genetic interaction data.
}
\description{
These functions implement a classifier to classify three classes of pathway membership of the RasMAPK and JNK pathway. For each sample and each channel a sparse linear discriminat classifier is trained. The posterior probabilities are averaged over all single classifiers. The classification posterior probabilities of three classes are plotted as a ternary plot (ternary plot adapted from CRAN package \code{vcd}).
}
\usage{
MAPK.cv.classifier(sgi, traingroups)
MAPK.predict.classification(sgi, traingroups)
MAPK.plot.classification(posterior,
                         classes = NULL, classnames = NULL, 
		         col = "darkgray", y = NULL,
                         classcol = NULL,
		         main = "predicted classification probabilities",
                         pop = TRUE, threshText = 0.3,
			 textToLeft = NULL, textToRight = NULL)
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}}
  \item{traingroups}{A list of gene names for the training examples. For each class there should be a vector of gene names.}
  \item{posterior}{A matrix of posterior probabilities. Each row represents one gene, each column represents one class.}
  \item{classes}{The three classes to be displayed on the ternary plot.}
  \item{classnames}{The class names to be displayed.}
  \item{col}{The color used for the text labels.}
  \item{y}{A factor representing the class label for each gene in posterior.}
  \item{classcol}{The color used for the three classes.}
  \item{main}{The title of the plot.}
  \item{pop}{If TRUE, all viewports are popped before finishing the function.}
  \item{threshText}{A threshold for the posterior probability of the three classes. Only genes that are assigned with a larger probability to the three classes are shown.}
  \item{textToLeft}{These text labels will be shown on the left hand side.}
  \item{textToRight}{These text labels will be shown on the right hand side.}
}
\details{
The code for the ternary plot (used by \code{MAPK.plot.classification}) is adapted from the function \code{ternaryplot} in the CRAN package \code{vcd}
Author of the original code is David Meyer (David.Meyer@R-project.org). References: M. Friendly (2000), Visualizing Categorical Data. SAS Institute, Cary, NC.
This code is specialized for the publication "Mapping Signalling Networks by RNAi ..." in Nat. Methods. It is highly recommended to use the original code by David Meyer.
}
\value{
\code{MAPK.cv.classifier} returns a list with the cross validated class assignment probability, as well as the results of the single classifiers.

\code{MAPK.predict.classifier} returns the predicted posterior probabilities of new genes as well as the classification results of the single classifier.

\code{MAPK.plot.classifier} returns nothing.
}
\references{
function \code{ternaryplot}, CRAN package \code{vdc}.
M. Friendly (2000), Visualizing Categorical Data. SAS Institute, Cary, NC.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract}}, \code{\link{RNAinteractMAPK}}
}
\keyword{ misc }
