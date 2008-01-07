\encoding{UTF-8}
\name{genind2genpop}
\alias{genind2genpop}
\title{Conversion from a genind to a genpop object}
\description{The function \code{genind2genpop} converts genotypes data
  (genind) into alleles counts per population (genpop).}
\usage{
genind2genpop(x,pop=NULL,missing=NA,quiet=FALSE)
}
\arguments{
  \item{x}{an object of class \code{genind}.}
  \item{pop}{a factor giving the population of each genotype in 'x'. If
    note provided, seeked in x@pop, but if given, the argument
    prevails on x@pop.}
  \item{missing}{can be NA, 0, or "replace". See details for more information.}
  \item{quiet}{logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
    
}
\details{The values of the 'missing' argument in \code{genind2genpop}
  have the following effects:\cr
  - NA: if all genotypes of a population for a given allele are
  missing, count value will be NA\cr
  
  - 0: if all genotypes of a population for a given allele are
  missing, count value will be 0\cr
  
  - "replace": when an allele is not typed in a population, it is
  assigned an allele count so that the allelic frequency in this
  populations is the same as the frequency in the whole dataset.\cr
  If allele 'j' of locus 'k' in pop 'i' is missing, the
  count value is number 'x' so that the frequency 'x/s' ('s' being the
  number of observations in 'k' ) equals the frequency 'f' computed on the whole data (i.e. considering all pop as one)\cr
  Then x verifies:\cr
  \deqn{x/s = f(1-f) => x=f(1-f)s}

  Note on 'pop': the factor is sorted in the order of appearance of each modality.
}
\value{
  A genpop object. The component @other in 'x' is passed to the created
  genpop object.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\seealso{
  \code{"\linkS4class{genind}"}, \code{"\linkS4class{genpop}"}
}
\examples{
  data(nancycats)
  nancycats
  catpop <- genind2genpop(nancycats)
  catpop
  summary(catpop)
}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}
