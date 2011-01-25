\encoding{UTF-8}
\name{read.snp}
\alias{read.snp}
\title{ Reading Single Nucleotide Polymorphism data}
\description{
  The function \code{read.snp} reads a SNP data file with extension '.snp' and
  converts it into a \linkS4class{genlight} object. This format is
  devoted to handle biallelic SNP only, but can accomodate massive
  datasets such as complete genomes with considerably less memory than
  other formats.

  A description of the .snp format is provided in an example file
  distributed with adegenet (see example below).
}
\usage{
read.snp(file, quiet=FALSE, \dots)
}
\arguments{
   \item{file}{ a character string giving the path to the file to
    convert, with the extension ".snp".}
  \item{quiet}{ logical stating whether a conversion messages should be
    printed (TRUE,default) or not (FALSE).}
  \item{\dots}{other arguments to be passed to other functions -
  currently not used.}
}
\details{
  Details of the .snp format can be found in the example file
  distributed with adegenet (see below), or on the adegenet website
  (type \code{adegenetWeb()} in R).
}
\value{an object of the class \linkS4class{genlight}}
\seealso{
\code{\link{import2genind}}, \code{\link{df2genind}}, \code{\link{read.genetix}}
  \code{\link{read.fstat}}, \code{\link{read.structure}}, \code{\link{read.genepop}}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## show the example file ##
## this is the path to the file:
system.file("files/exampleSnpDat.snp",package="adegenet")

## show its content:
file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))

## read the file
obj <- read.snp(system.file("files/exampleSnpDat.snp",package="adegenet"))
obj
as.matrix(obj)
ploidy(obj)
alleles(obj)
locNames(obj)
}
\keyword{manip}
