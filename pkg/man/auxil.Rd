\encoding{UTF-8}
\docType{methods}
\name{Auxiliary functions}
\alias{adegenetWeb}
\alias{adegenetTutorial}
\alias{checkType}
\alias{.rmspaces}
\alias{.genlab}
\alias{.readExt}
\alias{corner}
\alias{num2col}
\alias{transp}
\title{ Utilities functions for adegenet}
\description{
  These functions are mostly used internally in adegenet.  The notable
  exceptions are \code{adegenetWeb} which opens the adegenet website in
  the default navigator, \code{adegenetTutorial} which opens online
  tutorials for adegenet, and a few auxiliary graphical functions.\cr

  The other functions are:\cr
  - \code{checkType}: checks the type of markers being used in a
  function and issues an error if appropriate.\cr
  - \code{.rmspaces}: remove peripheric spaces in a character string. \cr
  - \code{.genlab}: generate labels in a correct alphanumeric ordering. \cr
  - \code{.readExt}: read the extension of a given file. \cr
  - \code{corner}: adds text to a corner of a figure. \cr
  - \code{num2col}: translates a numeric vector into colors. \cr
  - \code{transp}: adds transparency to a vector of colors. Note that
  transparent colors are not supported on some graphical devices.\cr
}
\usage{
adegenetWeb()
adegenetTutorial(which=c("general","spca"))
.genlab(base, n)
corner(text, posi="topleft",  inset=0.1, \dots)
num2col(x, col.pal=heat.colors, reverse=FALSE,
               x.min=min(x), x.max=max(x), na.col="green")
transp(col, alpha=.5)
}
\arguments{
  \item{which}{a character string being "general" or "spca", indicating
    which tutorial should be opened.}
  \item{base}{a character string forming the base of the labels}
  \item{n}{the number of labels to generate}
  \item{text}{a character string to be added to the plot}
  \item{posi}{a character matching any combinations of "top/bottom" and
    "left/right".}
  \item{inset}{a vector of two numeric values (recycled if needed)
    indicating the inset, as a fraction of the plotting region.}
  \item{\dots}{further arguments to be passed to \code{\link{text}}}
  \item{x}{a numeric vector}
  \item{col.pal}{a function generating colors according to a given palette.}
  \item{reverse}{a logical stating whether the palette should be
    inverted (TRUE), or not (FALSE, default).}
  \item{x.min}{the minimal value from which to start the color scale}
  \item{x.max}{the maximal value from which to start the color scale}
  \item{na.col}{the color to be used for missing values (NAs)}
  \item{col}{a vector of colors}
  \item{alpha}{a numeric value between 0 and 1 representing the alpha
  coefficient; 0: total transparency; 1: no transparency.}
}
\value{
  For \code{.genlab}, a character vector of size "n".
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## this opens the adegenet website
adegenetWeb()

## this opens the general tutorial for adegenet
adegenetTutorial()

## this opens the sPCA tutorial for adegenet
adegenetTutorial("spca")
}

.genlab("Locus-",11)
}
\keyword{manip}