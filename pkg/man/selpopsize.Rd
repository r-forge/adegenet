\encoding{UTF-8}
\docType{methods}
\name{selpopsize}
\alias{selpopsize}
\alias{selpopsize-methods}
\alias{selpopsize,ANY-method}
\alias{selpopsize,genind-method}
\title{ Select genotypes of well-represented populations}
\description{
  The function \code{selpopsize} checks the sample size of each population in
  a \linkS4class{genind} object and keeps only genotypes of populations
  having a given minimum size.
}
\usage{
\S4method{selpopsize}{genind}(x,pop=NULL,nMin=10)
}
\arguments{
  \item{x}{a \linkS4class{genind} object}
  \item{pop}{a vector of characters or a factor giving the population of
    each genotype in 'x'. If not provided, seeked from x\$pop.}
  \item{nMin}{the minimum sample size for a population to be
    retained. Samples sizes strictly less than \code{nMin} will be
    discarded, those equal to or greater than \code{nMin} are kept.}
}
\value{
  A  \linkS4class{genind} object.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{seploc}}, \code{\link{repool}}}
\examples{
data(microbov)

table(pop(microbov))
obj <- selpopsize(microbov, n=50)

obj
table(pop(obj))

}
\keyword{manip}