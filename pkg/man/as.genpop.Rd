\encoding{UTF-8}
\name{as.genpop}
\alias{as.genpop}
\alias{is.genpop}
\title{as.genpop and is.genpop}
\description{Direct conversion from a table of alleles counts as
  to a genpop object.}
\usage{
as.genpop(tab=NULL,prevcall=NULL)
is.genpop(x)
}
\arguments{
  \item{tab}{a pop x alleles matrix which terms are numbers of alleles,
    i.e. like in a genpop object}  
  \item{prevcall}{call of an object}
  \item{x}{an object}
}
\value{
  For {as.genpop}, a genpop object. For {is.genpop}, a logical.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\seealso{
  \code{"\linkS4class{genpop}"}, and \code{\link{genind2genpop}} for
  conversion from a genind to a genpop object.
}
\examples{
data(nancycats)
obj <- genind2genpop(nancycats)

# isolate one locus, fca77
obj <- seploc(obj,truenames=TRUE)$"fca77"
class(obj)
dim(obj)
head(obj)
# make it a new genpop
obj <- as.genpop(obj)
obj
}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}
