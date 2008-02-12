\name{repool}
\alias{repool}
\title{Pool several genotypes into the same dataset}
\description{
  The function \code{repool} allows to merge genotypes in different
  \linkS4class{genind} objects into a single  \linkS4class{genind}.
  The argument of the function is a list of \linkS4class{genind} objects.
  The markers have to be the same for all objects to be merged.
  The alleles can differe.
}
\usage{
repool(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a list whose components are valid \linkS4class{genind} objects}
}
\value{
  A \linkS4class{genind} object.
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{seploc}}, \code{\link{seppop}}}
\examples{
}
\keyword{manip}
