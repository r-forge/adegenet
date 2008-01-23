\encoding{UTF-8}
\name{import}
\alias{df2genind}
\alias{import2genind}
\alias{read.genetix}
\alias{read.genepop}
\alias{read.fstat}
\title{ Conversion functions for adegenet}
\description{
  Their are two ways to import genotyping data to a \linkS4class{genind} object: from
  a data.frame with a given format (see below), or from a file with a
  recognized extension.\cr
  
  The function \code{df2genind} converts a data.frame (or a matrix) into
  a \linkS4class{genind} object. The data.frame must meet the following
  requirements:\cr
  - genotypes are in row (on row per genotype)\cr
  - markers are in columns\cr
  - each element is a string of (2, 4 or 6) characters coding both alleles without
  separation. Missing values are coded by NA or zeros. Uncomplete
  strings are filled with "0" at the begining. For instance, "912" in a
  4-characters coding scheme is interpreted as "0912", thus as an
  heterozygote with alleles "09" and "12". Note that this format is the
  one of GENETIX for a 6-characters coding scheme\cr
  
  The function \code{import2genind} detects the extension of
  the file given in argument and seeks for an appropriate import
  function to create a \code{genind} object.\cr
  Current functions are :\cr
  - \code{read.genetix} for GENETIX files (.gtx). Note that this
  function is called by the others. \cr
  - \code{read.genepop} for Genepop files (.gen) \cr
  - \code{read.fstat} for Fstat files \code{.dat} \cr
}
\usage{
df2genind(X, ncode=NULL, ind.names=NULL, pop=NULL, missing=NA)
import2genind(file,missing=NA,quiet=FALSE)
read.genetix(file=NULL,missing=NA,quiet=FALSE)
read.genepop(file,missing=NA,quiet=FALSE)
read.fstat(file,missing=NA,quiet=FALSE)
}
\arguments{
  \item{X}{a matrix or a data.frame (see decription)}
  \item{ncode}{an optional integer the number of characters used for coding one
    genotype at one locus (can be 2, 4 or 6). If not provided, this is
    determined from data.}
  \item{ind.names}{an optional character vector giving the individuals names}
  \item{pop}{an optional factor giving the population of each individual.}
  \item{file}{ a character string giving the path to the file to
    convert, with the appropriate extension.}
  \item{missing}{can be NA, 0 or "mean". See details section.}
  \item{quiet}{ logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
}
   
\details{There are 3 treatments for missing values: \cr
  - NA: kept as NA.\cr
  
  - 0: allelic frequencies are set to 0  on all alleles of the concerned
  locus. Recommended for a PCA on compositionnal data.\cr
  
  - "mean": missing values are replaced by the mean frequency of the
  corresponding allele, computed on the whole set of
  individuals. Recommended for a centred PCA.\cr

  Beware: same data in different formats are not expected to produce 
  exactly the same \code{genind} objects.\cr
  For instance, conversions made by GENETIX to Fstat may change the
  the sorting of the genotypes; GENETIX stores individual names whereas
  Fstat does not; Genepop chooses a sample's name from the name of its last
  genotype; etc.
}
\value{an object of the class \code{genind}}
\references{
  Belkhir K., Borsa P., Chikhi L., Raufaste N. & Bonhomme F. (1996-2004)
  GENETIX 4.05, logiciel sous Windows TM pour la génétique des
  populations. Laboratoire Génome, Populations, Interactions, CNRS UMR
  5000, Université de Montpellier II, Montpellier (France). \cr
  
  Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population
  genetics software for exact tests and ecumenicism. \emph{J. Heredity},
  \bold{86}:248-249 \cr
  
  Fstat (version 2.9.3). Software by Jerome
  Goudet. http://www2.unil.ch/popgen/softwares/fstat.htm\cr

  Excoffier L. & Heckel G.(2006) Computer programs for population genetics data analysis: a survival guide \emph{Nature}, \bold{7}: 745-758
}
\seealso{\code{\link[pkg:hierfstat]{read.fstat.data}} in package hierfstat}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
read.genetix(system.file("files/nancycats.gtx",package="adegenet"))

read.fstat(system.file("files/nancycats.dat",package="adegenet"))

read.genepop(system.file("files/nancycats.gen",package="adegenet"))

import2genind(system.file("files/nancycats.gtx",
package="adegenet"))

if(require(hierfstat)){
obj <- read.fstat(system.file("data/diploid.dat",package="hierfstat"))
obj
}
}
\keyword{manip}
