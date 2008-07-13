\name{scaleGen-methods}
\docType{methods}
\alias{scaleGen}
\alias{scaleGen-methods}
\alias{scaleGen,genind-method}
\alias{scaleGen,genpop-method}
\title{ Compute scaled allele frequencies }
\description{
  The generic function \code{scaleGen} is an analogue to the
  \code{scale} function, but is designed with further arguments giving
  scaling options.\cr

  Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
  objects.
  Both return data.frames of scaled allele frequencies.
}
\usage{
\S4method{scaleGen}{genind}(x, center=TRUE, scale=TRUE, method=c("sigma", "binom"), missing=c("NA","0","mean"),truenames=TRUE)
\S4method{scaleGen}{genpop}(x, center=TRUE, scale=TRUE, method=c("sigma", "binom"), missing=c("NA","0","mean"),truenames=TRUE)
}
\arguments{
  \item{x}{a \linkS4class{genind} and \linkS4class{genpop} object}
  \item{center}{a logical stating whether alleles frequencies should be
    centred to mean zero (default to TRUE). Alternatively, a vector of
    numeric values, one per allele, can be supplied: these values will be
    substracted from the allele frequencies.}
  \item{scale}{a logical stating whether alleles frequencies should be
    scaled (default to TRUE). Alternatively, a vector of
    numeric values, one per allele, can be supplied: these values will be
    substracted from the allele frequencies.}
  \item{method}{a character indicating the method to be used. See details.}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
 \item{missing}{a character giving the treatment for missing values. Can be "NA", "0" or "mean"}
 }
 \value{
   A \linkS4class{genind} and \linkS4class{genpop} object.
}
\details{
  The argument \code{method} is used as follows:\cr

  - \code{sigma}: scaling is made using the usual standard deviation\cr

  - \code{binom}: scaling is made using the theoretical variance of the
  allele frequency. This can be used to avoid that frequencies close to
  0.5 have a stronger variance that those close to 0 or 1.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
## load data
data(microbov)
obj <- genind2genpop(microbov)

## compare different scaling
X1 <- scaleGen(obj)
X2 <- scaleGen(obj,met="bin")

if(require(ade4)){
## compute PCAs
pcaObj <- dudi.pca(obj,scale=FALSE,scannf=FALSE) # pca with no scaling
pcaX1 <- dudi.pca(X1,scale=FALSE,scannf=FALSE) # pca with usual scaling
pcaX2 <- dudi.pca(X2,scale=FALSE,scannf=FALSE) # pca with scaling for binomial variance

## change signs of loadings if needed
if(cor(pcaX1$c1[,1],pcaX2$c1[,1]) < 0) {pcaX1$c1[,1] <- -pcaX1$c1[,1]}
if(cor(pcaX1$c1[,2],pcaX2$c1[,2]) < 0) {pcaX1$c1[,2] <- -pcaX1$c1[,2]}

## graphics
par(mfrow=c(2,2))
barplot(pcaObj$eig,main="Eigenvalues\n no scaling")
barplot(pcaX1$eig,main="Eigenvalues\n usual scaling")
barplot(pcaX2$eig,main="Eigenvalues\n 'binomial' scaling")
s.match(pcaX1$c1,pcaX2$c1,clab=0,sub="usual -> binom")
}
}
\keyword{methods}
\keyword{manip}
