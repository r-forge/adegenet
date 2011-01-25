\encoding{UTF-8}
\name{genlight auxiliary functions}
\alias{glSum}
\alias{glNA}
\alias{glMean}
\alias{glVar}
\title{Auxiliary functions for genlight objects}
\description{
  These functions provide facilities for usual computations using
  \linkS4class{genlight} objects. In many cases, the output
  depends on whether the information units are individuals, or
  alleles within individuals (see details).
  
  These functions are:
  
  - \code{glSum}: computes the sum of the number of second allele in each SNP.

  - \code{glNA}: computes the number of missing values in each SNP.

  - \code{glMean}: computes the mean number of second allele in each SNP.

  - \code{glVar}: computes the variance of the number of second allele in each SNP.
}
\usage{
glSum(x, alleleAsUnit = TRUE)
glNA(x, alleleAsUnit = TRUE)
glMean(x, alleleAsUnit = TRUE)
glVar(x, alleleAsUnit = TRUE)
}
\arguments{
  \item{x}{a \linkS4class{genlight} object}
  \item{alleleAsUnit}{a logical indicating whether alleles are
    considered as units (i.e., a diploid genotype equals two samples, a
    triploid, three, etc.) or whether individuals are considered as
    units of information.}
}
\details{
  === On the unit of information ===
  In the cases where individuals can have different ploidy, computation
  of sums, means, etc. of allelic data depends on what we consider as a
  unit of information.

  To estimate e.g. allele frequencies, unit of information can be
  considered as the allele, so that a diploid genotype contains two
  samples, a triploid individual, three samples, etc. In such a case,
  all computations are done directly on the number of alleles. This
  corresponds to \code{alleleAsUnit = TRUE}.

  However, when the focus is put on studying differences/similarities
  between individuals, the unit of information is the individual, and
  all genotypes possess the same information no matter what their ploidy
  is. In this case, computations are made after standardizing
  individual genotypes to relative allele frequencies. This
  corresponds to \code{alleleAsUnit = FALSE}.

  Note that when all individuals have the same ploidy, this distinction
  does not hold any more.
}
\value{
 A numeric vector containing the requested information.
}

\seealso{
  - \code{\link{genlight}}: class of object for storing massive binary
  SNP data.
  
  - \code{\link{dapc}}: Discriminant Analysis of Principal Components.

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(2,1,1,1,1,NA)))
x
as.matrix(x)
ploidy(x)

## compute statistics - allele as unit ##
glNA(x)
glSum(x)
glMean(x)

## compute statistics - individual as unit ##
glNA(x, FALSE)
glSum(x, FALSE)
glMean(x, FALSE)

## explanation: data are taken as relative frequencies
temp <- as.matrix(x)/ploidy(x)
apply(temp,2, function(e) sum(is.na(e))) # NAs
apply(temp,2,sum, na.rm=TRUE) # sum
apply(temp,2,mean, na.rm=TRUE) # mean

}
\keyword{multivariate}