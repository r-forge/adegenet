\name{as methods in adegenet}
\docType{methods}
\alias{as-method}
\alias{as,genind,data.frame-method}
\alias{as,genpop,data.frame-method}
\alias{as,genind,matrix-method}
\alias{as,genpop,matrix-method}
\alias{as,genind,genpop-method}
\alias{coerce,genind,data.frame-method}
\alias{coerce,genpop,data.frame-method}
\alias{coerce,genind,matrix-method}
\alias{coerce,genpop,matrix-method}
\alias{coerce,genind,genpop-method}
\alias{as.data.frame.genind}
\alias{as.data.frame.genpop}
\alias{as.matrix.genind}
\alias{as.matrix.genpop}
\alias{as.genpop.genind}
\title{Converting genind/genpop objects to other classes}
\section{Usage}{
\code{as(object, Class)}
}
\section{Arguments}{
\describe{
  \item{\code{object}}{a \linkS4class{genind} or a \linkS4class{genpop} object.}
  \item{\code{Class}}{the name of the class to which the object should
    be coerced, for instance \code{"data.frame"} or \code{"matrix"}.}
}
}
\description{
 These S3 and S4 methods are used to coerce \linkS4class{genind} and
 \linkS4class{genpop} objects to matrix-like objects. In most cases,
 this is equivalent to calling the \code{@tab} slot.
}
\section{Methods}{
  \item{coerce}{from one object class to another using
    \code{as(object,"Class")}, where the \code{object} is of the old
    class and the returned object is of the new class \code{"Class"}.
  }
}

\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\examples{
data(microbov)
x <- na.replace(microbov,method="0")
as(x[1:3],"data.frame")
}
\keyword{methods}