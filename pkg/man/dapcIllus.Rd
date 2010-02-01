\encoding{UTF-8}
\name{dapcIllus}
\alias{dapcIllus}
\docType{data}
\title{Simulated data illustrating the dapc}
\description{
  Datasets illustrating the Discriminant Analysis of Principal
  Components (DAPC, Jombart et al. submitted).\cr

  These data were simulated using various models using Easypop (2.0.1).
  The \code{dapcIllus} is a list containing the
  following \linkS4class{genind} objects:\cr
  - "a": island model with 6 populations \cr
  - "b": hierarchical island model with 6 populations (3,2,1) \cr
  - "c": one-dimensional stepping stone with 2x6 populations, and a boundary between the
  two sets of 6 populations\cr
  - "d": one-dimensional stepping stone with 24 populations\cr

  See "source" for a reference providing simulation details.
}
\usage{data(dapcIllus)}
\format{
  \code{dapcIllus} is list of 4 components being all genind objects.
}
\source{
  Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{PLoS genetics}.

}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\seealso{ \code{\link{dapc}}, \code{\link{eHGDP}},  \code{\link{H3N2}}}
\references{
   Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{PLoS genetics}.
}
\examples{
if(require(MASS) & require(ade4)){

data(dapcIllus)
attach(dapcIllus)

## PERFORM DAPCs
dapc.a <- dapc(dapcIllus$a, n.pca=100, n.da=5)
dapc.b <- dapc(dapcIllus$b, n.pca=100, n.da=5)
dapc.c <- dapc(dapcIllus$c, n.pca=100, n.da=11)
dapc.d <- dapc(dapcIllus$d, n.pca=100, n.da=11)


## LOOK AT ONE RESULT
dapc.a
summary(dapc.a)

## FORM A LIST OF RESULTS FOR THE 4 DATASETS
lres <- list(dapc.a, dapc.b, dapc.c, dapc.d)


## DRAW 4 SCATTERPLOTS
par(mfrow=c(2,2))
sapply(lres, scatter)


# detach data
detach(dapcIllus)
}
}
\keyword{datasets}
