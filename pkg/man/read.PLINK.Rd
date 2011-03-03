\encoding{UTF-8}
\name{read.PLINK}
\alias{read.PLINK}
\title{ Reading PLINK Single Nucleotide Polymorphism data}
\description{
  The function \code{read.PLINK} reads a data file exported by the PLINK
  software with extension '.raw' and converts it into a
  \linkS4class{genlight} object.

  The function reads data by chunks of several genomes (minimum 1, no
  maximum) at a time, which allows one to read massive datasets with
  negligible RAM requirements (albeit at a cost of computational
  time). The argument \code{chunkSize} indicates the number of genomes
  read at a time. Increasing this value decreases the computational time
  required to read data in, while increasing memory requirements.

  See details for the documentation about how to export data using PLINK
  to the '.raw' format.
}
\usage{
read.PLINK(file, quiet=FALSE, chunkSize=1000, ploidy=2,
           multicore=require("multicore"), n.cores=NULL, \dots)
}
\arguments{
  \item{file}{ a character string giving the path to the file to
    convert, with the extension ".PLINK".}
  \item{quiet}{ logical stating whether a conversion messages should be
    printed (TRUE,default) or not (FALSE).}
  \item{chunkSize}{an integer indicating the number of genomes to be
    read at a time; larger values require more RAM but decrease the time
    needed to read the data.}
  \item{ploidy}{an integer indicating the ploidy of the data; defaults to
    2.}
  \item{multicore}{a logical indicating whether multiple cores -if
    available- should be used for the computations (TRUE, default), or
    not (FALSE); requires the package \code{multicore} to be installed
    (see details).}
  \item{n.cores}{if \code{multicore} is TRUE, the number of cores to be
    used in the computations; if NULL, then the maximum number of cores
    available on the computer is used.}
  \item{\dots}{other arguments to be passed to other functions -
    currently not used.}
}
\details{
  === Exporting data from PLINK ===
  
  Data need to be exported from PLINK using the option "--recodeA" (and
  NOT "--recodeAD"). The PLINK command should therefore look like:
  \code{plink --file data --recodeA}. For more information on this
  topic, please look at this webpage:
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml}

  
  === Using multiple cores ===
  
  Most recent machines have one or several processors with multiple
  cores. R processes usually use one single core. The package
  \code{multicore} allows for parallelizing some computations on
  multiple cores, which decreases drastically computational time.

  To use this functionality, you need to have the last version of the
  \code{multicore} package installed. To install it, type:
  install.packages("multicore",,"http://rforge.net/",type="source")

  DO NOT use the version on CRAN, which is slightly outdated.
}
\value{an object of the class \linkS4class{genlight}}
\seealso{
  - a description of the class \linkS4class{genlight}.

  - \code{\link{read.snp}}: read from adegenet's dedicated format to
  - store biallelic SNP data.
  
  - other import function in adegenet: \code{\link{import2genind}},
  \code{\link{df2genind}}, \code{\link{read.genetix}}
  \code{\link{read.fstat}}, \code{\link{read.structure}},
  \code{\link{read.genepop}}.

  - another function \code{read.plink} is available in the package
  \code{snpMatrix}.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{

}
\keyword{manip}
