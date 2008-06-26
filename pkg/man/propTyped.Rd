\name{propTyped-methods}
\docType{methods}
\alias{propTyped}
\alias{propTyped-methods}
\alias{propTyped,genind-method}
\alias{propTyped,genpop-method}
\title{ Compute the proportion of typed elements }
\description{
  The generic function \code{propTyped} is devoted to investigating the
  structure of missing data in adegenet objects.\cr

  Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
  objects. They can return the proportion of available
  (i.e. non-missing) data per individual/population, locus, or the
  combination of both in with case the matrix indicates which entity
  (individual or population) was typed on which locus.
}
\usage{
\S4method{propTyped}{genind}(x,  by=c("ind","loc","both"))
\S4method{propTyped}{genpop}(x,  by=c("pop","loc","both"))
}
\arguments{
  \item{x}{a \linkS4class{genind} and \linkS4class{genpop} object}
  \item{by}{a character being "ind","loc", or "both" for
  \linkS4class{genind} object and "pop","loc", or "both" for
  \linkS4class{genpop} object. It specifies whether proportion of typed
  data are provided by entity ("ind"/"pop"), by locus ("loc") or both
  ("both"). See details.}
 }
 \value{
   A vector of proportion (when \code{by} equals "ind", "pop", or
  "loc"), or a matrix of binary data (when \code{by} equals "both")
}
\details{
  The argument \code{method} is used as follows:\cr

  - \code{sigma}: scaling is made using the usual standard deviation\cr

  - \code{binom}: scaling is made using the theoretical variance of the
  allele frequency. This can be used to avoid that frequencies close to
  0.5 have a stronger variance that those close to 0 or 1.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{

}
\keyword{methods}
\keyword{manip}
