\encoding{UTF-8}
\docType{methods}
\name{Auxiliary functions}
\alias{truenames}
\alias{truenames-methods}
\alias{truenames,ANY-method}
\alias{truenames,genind-method}
\alias{truenames,genpop-method}
\alias{seploc}
\alias{seploc-methods}
\alias{seploc,ANY-method}
\alias{seploc,genind-method}
\alias{seploc,genpop-method}
\alias{.rmspaces}
\alias{.readExt}
\alias{adegenetWeb}
\title{ Utilities functions}
\description{ Most of theses functions are S4 methods to be used with \code{genind} and
  \code{genpop} objects.\cr
   The function \code{seploc} splits an object by marker, returning a
   list whose components each correspond to a marker.\cr
   
   \code{truenames} returns elements of the object using true names
  (as opposed to generic labels) for individuals, markers, alleles, and population.\cr
  An exception is \code{adegenetWeb}, which simply opens the adegenet
  website in the default navigator.
}
\usage{
\S4method{seploc}{genind}(truenames=FALSE,res.type=c("genind","matrix"))
\S4method{seploc}{genpop}(truenames=FALSE,res.type=c("genpop","matrix"))
adegenetWeb()
}
\arguments{
    \item{truenames}{a logical indicating whether generic labels
      (FALSE,default) or true names should be used (TRUE).}
    \item{res.type}{a character indicating the type of returned results,
      a genind or genpop object (default) or a matrix of data.}
}
\value{The function \code{seploc} returns an list of objects of the
  same class as the initial object, or a list of matrices similar to
  x\$tab.\cr

  \code{truenames} returns a matrix similar to the x\$tab slot but with
  true labels. If x\$pop exists, it returns a list with this
  matrix (\$tab) and a population vector with true names (\$pop).\cr
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
data(microbov)
# restore true names
truenames(microbov)$tab[1:5,1:5]

# separate all markers
obj <- seploc(microbov)
names(obj)

obj$INRA5
}
\keyword{manip}