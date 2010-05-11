\name{haploGen}
\alias{haploGen}
\alias{print.haploGen}
\alias{[.haploGen}
\alias{labels.haploGen}
\alias{as.POSIXct.haploGen}
\alias{seqTrack.haploGen}
\alias{haploGen-class}
% \alias{seqTrackG.haploGen}
% \alias{optimize.seqTrack.haploGen}
\alias{as.seqTrack.haploGen}
\alias{plotHaploGen}
\alias{sample.haploGen}
\alias{as,haploGen,graphNEL-method}
\alias{coerce,haploGen,graphNEL-method}
\title{Simulation of populations of haplotypes}
\description{
  Important: these functions are parts of a publication currently under
  review. They will be documented once accepted for publication.
  Please email the author if you are interested in using it.
}
\usage{
haploGen(seq.length=1000, mu=0.0001, t.max=50,
                     gen.time=function(){round(rnorm(1,5,1))},
                     repro=function(){round(rnorm(1,2,1))}, max.nb.haplo=1e3,
                     geo.sim=TRUE, grid.size=5, lambda.xy=0.5,
                     mat.connect=NULL,
                     ini.n=1, ini.xy=NULL)

print.haploGen(x, ...)

"[.haploGen"(x,i,j,drop=FALSE)

labels.haploGen(object, ...)

as.POSIXct.haploGen(x, tz="", origin=as.POSIXct("2000/01/01"), ...)

seqTrack.haploGen(x, best=c("min","max"), prox.mat=NULL, ...)

as.seqTrack.haploGen(x)

plotHaploGen(x, annot=FALSE, date.range=NULL, col=NULL, bg="grey", add=FALSE, ...)

sample.haploGen(x, n)

setAs("haploGen", "graphNEL", def=function(from)

}
\arguments{
  \item{seq.length}{}
  \item{mu}{}
  \item{t.max}{}
  \item{gen.time}{}
  \item{repro}{}
  \item{max.nb.haplo}{}
  \item{geo.sim}{}
  \item{grid.size}{}
  \item{lambda.xy}{}
  \item{mat.connect}{}
  \item{ini.n}{}
  \item{ini.xy}{}
  \item{x,object}{}
  \item{i,j}{}
  \item{best}{}
  \item{prox.mat}{}
  \item{annot,date.range,col,bg,add}{}
  \item{n}{}
 }
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
