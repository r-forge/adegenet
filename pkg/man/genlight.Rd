\name{genlight-class}
\docType{class}
\alias{genlight}
\alias{genlight-class}
\alias{[,genlight-method}
\alias{[,genlight,ANY,ANY-method}
\alias{initialize,genlight-method}
\alias{show,genlight-method}
\alias{nLoc,genlight-method}
\alias{nInd,genlight-method}
\alias{$,genlight-method}
\alias{$<-,genlight-method}
\alias{names,genlight-method}
\alias{ploidy,genlight-method}
\alias{locNames,genlight-method}
\alias{indNames,genlight-method}
\alias{alleles,genlight-method}
\alias{pop,genlight-method}
\alias{pop<-,genlight-method}
\alias{as,genlight,matrix-method}
\alias{coerce,genlight,matrix-method}
\alias{as.matrix.genlight}
\alias{as,genlight,data.frame-method}
\alias{coerce,genlight,data.frame-method}
\alias{as.data.frame.genlight}
\alias{as,genlight,list-method}
\alias{coerce,genlight,list-method}
\alias{as.list.genlight}
% \alias{,genlight-method}
% \alias{,genlight-method}
% \alias{,genlight-method}
% \alias{,genlight-method}
%%%%
\title{Formal class "genlight"}
\description{
  The class \code{genlight} is a formal (S4) class for storing a genotypes
  of binary SNPs in a compact way, using a bit-level coding scheme.
  This storage is most efficient with haploid data, where the memory
  taken to represent data can reduced more than 50 times. However,
  \code{genlight} can be used for any level of ploidy, and still remain an
  efficient storage mode.

  A \code{genlight} object can be constructed from vectors of integers
  giving the number of the second allele for each locus and each
  individual (see 'Objects of the class genlight' below).

  \code{genlight} stores a multiple genotypes. Each genotype is stored
  as a \linkS4class{SNPbin} object.
}
\section{Objects from the class genlight}{
  \code{genlight} objects can be created by calls to \code{new("genlight",
    ...)}, where '...' can be the following arguments:
  
  \describe{
    \item{\code{gen}}{input genotypes, where each genotype is coded as a
      vector of numbers of the second allele. If a list, each slot of the
      list correspond to an individual; if a matrix or a data.frame, rows
      correspond to individuals and columns to SNPs. If individuals or
      loci are named in the input, these names will we stored in the
      produced object. All individuals are expected to have the same
      number of SNPs. Shorter genotypes are completed with NAs, issuing a
      warning.}
    \item{\code{ploidy}}{an optional vector of integers indicating the ploidy of the
      genotypes. Genotypes can therefore have different ploidy. If not
      provided, ploidy will be guessed from the data (as the
      maximum number of second alleles in each individual).}
    \item{\code{ind.names}}{an optional vector of characters giving the labels
      of the genotypes.}
    \item{\code{loc.names}}{an optional vector of characters giving the labels
      of the SNPs.}
    \item{\code{loc.all}}{an optional vector of characters indicating
    the alleles of each SNP; for each SNP, alleles must be coded by two
    letters separated by '/', e.g. 'a/t' is valid, but 'a  t' or 'a |t' are not.}
  }
}
\section{Slots}{
  The following slots are the content of instances of the class
  \code{genlight}; note that in most cases, it is better to retrieve
  information via accessors (see below), rather than by accessing the
  slots manually.
  \describe{
    \item{\code{gen}:}{a list of genotypes stored as  \linkS4class{SNPbin} objects.}
    \item{\code{n.loc}:}{an integer indicating the number of SNPs of the
      genotype.}
    \item{\code{ind.names}:}{a vector of characters indicating the names of
      genotypes.}
    \item{\code{loc.names}:}{a vector of characters indicating the names of
      SNPs.}
    \item{\code{loc.all}:}{a vector of characters indicating the alleles
      of each SNP.}
    \item{\code{ploidy}:}{a vector of integers indicating the ploidy of each genotype.}
  }
}
\section{Methods}{
  Here is a list of methods available for \code{genlight} objects. Most of
    these methods are accessors, that is, functions which are used to
    retrieve the content of the object. Specific manpages can exist for
    accessors with more than one argument. These are indicated by a '*'
    symbol next to the method's name. This list also contains methods
    for conversion from \code{genlight} to other classes.
  \describe{
    \item{[}{\code{signature(x = "genlight")}: usual method to subset
      objects in R. Is to be applied as if the object was a matrix where
      genotypes are rows and SNPs are columns. Indexing can be done via
      vectors of signed integers or of logicals.}
    \item{show}{\code{signature(x = "genlight")}: printing of the
      object.}
    \item{$}{\code{signature(x = "genlight")}: similar to the @ operator;
      used to access the content of slots of the object.}
    \item{$<-}{\code{signature(x = "genlight")}: similar to the @ operator;
      used to replace the content of slots of the object.}
    \item{nInd}{\code{signature(x = "genlight")}: returns the number of
      individuals in the object.}
    \item{nLoc}{\code{signature(x = "genlight")}: returns the number of
      SNPs in the object.}
    \item{names}{\code{signature(x = "genlight")}: returns the names of
      the slots of the object.}
    \item{ploidy}{\code{signature(x = "genlight")}: returns the ploidy of
      the genotypes.}
    \item{indNames}{\code{signature(x = "genlight")}: returns the names of
      the individuals, if provided when the object was contructed.}
    \item{locNames}{\code{signature(x = "genlight")}: returns the names of
      the loci, if provided when the object was contructed.}
    \item{alleles}{\code{signature(x = "genlight")}: returns the names
      of the alleles of each SNPs, if provided when the object was
      contructed.}
    \item{pop}{\code{signature(x = "genlight")}: returns a factor
      indicating the population of each individual, if provided when the
      object was contructed.}
    \item{pop<-}{\code{signature(x = "genlight")}: sets the population
    of each individual using a factor of length \code{nInd(x)}.}  
    \item{as.matrix}{\code{signature(x = "genlight")}: converts a
      \code{genlight} object into a matrix of integers, with individuals
      in rows and SNPs in columns. The S4 method 'as' can be used as
      well (e.g. as(x, "matrix")).}
    \item{as.data.frame}{\code{signature(x = "genlight")}: same as \code{as.matrix}.}
    \item{as.list}{\code{signature(x = "genlight")}: converts a
      \code{genlight} object into a list of genotypes coded as vector of
      integers (numbers of second allele). The S4 method 'as' can be
      used as well (e.g. as(x, "list")).}
  }
}
\author{Thibaut Jombart (\email{t.jombart@imperial.ac.uk})}
\seealso{
 Related class:\cr
  -  \code{\linkS4class{SNPbin}}, for storing individual genotypes of
  binary SNPs\cr
  
  -  \code{\linkS4class{genind}}, for storing other types of genetic markers. \cr
}
\examples{
## TOY EXAMPLE ##
## create and convert data
dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
x <- new("genlight", dat)
x

## examine the content of the object
names(x)
x@gen
x@gen[[1]]@snp # bit-level coding for first individual

## conversions
as.list(x)
as.matrix(x)

## round trips - must return TRUE
identical(x, new("genlight", as.list(x))) # list
identical(x, new("genlight", as.matrix(x))) # matrix
identical(x, new("genlight", as.data.frame(x))) # data.frame

## test subsetting
x[c(1,3)] # keep individuals 1 and 3
as.list(x[c(1,3)])
x[c(1,3), 1:2] # keep individuals 1 and 3, loci 1 and 2
as.list(x[c(1,3), 1:2])
x[c(TRUE,FALSE), c(TRUE,TRUE,FALSE,FALSE)] # same, using logicals
as.list(x[c(TRUE,FALSE), c(TRUE,TRUE,FALSE,FALSE)])


## REAL-SIZE EXAMPLE ##
## 50 genotypes of 1,000,000 SNPs
dat <- lapply(1:50, function(i) sample(c(0,1,NA), 1e6, prob=c(.5, .49, .01), replace=TRUE))
names(dat) <- paste("indiv", 1:length(dat))
print(object.size(dat), unit="aut") # size of the original data

x <- new("genlight", dat) # conversion
x
print(object.size(x), unit="au") # size of the genlight object
object.size(dat)/object.size(x) # conversion efficiency
}
\keyword{classes}
