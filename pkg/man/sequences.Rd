\encoding{UTF-8}
\name{SequencesToGenind}
\alias{DNAbin2genind}
\title{ Importing data from an alignement of sequences to a genind object}
\description{
  These functions take an alignement of sequences and translate SNPs
  into a \linkS4class{genind} object. Note that only polymorphic loci
  are retained.\cr

  Currently, accepted sequence formats are:\cr
  - DNAbin (ape package): function DNAbin2genind\cr
  - alignement (seqinr package): to come...
}
\usage{
DNAbin2genind(x, pop=NULL, na.char=c("n","-","?"))
}
\arguments{
 \item{x}{an object containing aligned sequences.}
 \item{pop}{an optional factor giving the population to which each sequence belongs.}
 \item{na.char}{a character vector providing values that should be
 considered as NA.}
}

\value{an object of the class \linkS4class{genind}}

\seealso{\code{\link{import2genind}}, \code{\link{read.genetix}},
  \code{\link{read.fstat}}, \code{\link{read.structure}},
  \code{\link{read.genepop}}, \code{\link[pkg:ape]{DNAbin}}.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
if(require(ape)){
data(woodmouse)
x <- DNAbin2genind(woodmouse)
x
genind2df(x)
}
}
\keyword{manip}
