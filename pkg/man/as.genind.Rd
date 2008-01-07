\encoding{UTF-8}
\name{as.genind}
\alias{as.genind}
\alias{is.genind}
\title{as.genind and is.genind}
\description{Auxiliary functions for genind objects.\cr
  \code{as.genind} performs a direct conversion from a table of
  genotypes coded as allelic frequency table to a genind object.
  \code{is.genind} tests if an object is a valid genind object.
}

\usage{
as.genind(tab=NULL,pop=NULL,prevcall=NULL)
is.genind(x)
}

\arguments{
  \item{tab}{a individuals x alleles matrix of genotypes coded as
    allelic frequencies, i.e. like in a genind object}
  \item{pop}{a factor giving the population of each genotype in 'x'}
  \item{prevcall}{call of an object}
  \item{x}{an object}
}
\value{
  For {as.genind}, a genind object. For {is.genind}, a logical.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\seealso{
  \code{"\linkS4class{genind}"}, and \code{\link{import2genind}} for
  import from various types of file.
}
\examples{
data(nancycats)
nancycats@loc.names

# isolate one marker, fca23
obj <- seploc(nancycats,truenames=TRUE)$"fca23"
class(obj)
dim(obj)

# make it a new genind
obj <- as.genind(obj)
obj
}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}
