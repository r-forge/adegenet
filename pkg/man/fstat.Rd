\encoding{UTF-8}
\name{fstat}
\alias{fstat}
\alias{Fst}
\alias{FST}
\alias{fst}
\alias{pairwise.fst}
\title{F statistics for genind objects}
\description{
  The function \code{fstat} computes a global Fst, while
  \code{pairwise.fst} computes Nei's pairwise Fst between all pairs of
  populations. Both functions are designed for \linkS4class{genind}
  objects.

  \code{fstat} is wrapper for \code{varcomp.glob} from package
  \code{hierfstat} for \linkS4class{genind} objects. It computes F statistics (Fst, Fis, Fit)
  given a set of genotypes and a grouping factor.

  \code{pairwise.fst} is an implementation of Nei's Fst in which
  heretozygosities are weighted by group sizes (see details).
}
\usage{
fstat(x, pop=NULL, fstonly=FALSE)

pairwise.fst(x, pop=NULL, res.type=c("dist","matrix"), truenames=TRUE)
}
\arguments{
  \item{x}{an object of class \linkS4class{genind}.}
  \item{pop}{a factor giving the 'population' of each individual. If NULL,
    pop is seeked from x\@pop. Note that the term population refers in
    fact to any grouping of individuals'.}
  \item{fstonly}{a logical stating whether only the Fst value should be
    returned (TRUE) instead of all F statistics (FALSE, default).}
  \item{res.type}{the type of result required: a \code{dist} object, or a
    symmetric matrix}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
}
\value{A matrix of F statistics.}
\seealso{\code{\link{Hs}}, \code{\link[hierfstat]{varcomp.glob}},
  \code{\link{gstat.randtest}}
}
\references{
  Nei, M. (1973) Analysis of gene diversity in subdivided
  populations. Proc Natl Acad Sci USA, 70: 3321-3323
}
\details{
  Let A and B be two populations of population sizes \eqn{n_A} and
  \eqn{n_B}, with expected heterozygosity Hs(A) and Hs(B),
  respectively. We denote Ht the expected heterozygosity of a population
  pooling A and B. Then, the pairwise Fst \eqn{Fst(A,B)} is computed
  as:\cr

  \eqn{\frac{(Ht - (n_A Hs(A) + n_B Hs(B))/(n_A + n_B) )}{Ht}} \cr
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
if(require(hierfstat)){
data(nancycats)

## Fst, Fis, Fit
fstat(nancycats)

## pairwise Fst
mat.fst <- pairwise.fst(nancycats, res.type="matrix")
mat.fst
}
}
\keyword{multivariate}