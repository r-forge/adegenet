\encoding{UTF-8}
\name{hybridize}
\alias{hybridize}
\title{ Simulated hybridization between two samples of populations}
\description{
  The function \code{hybridize} performs hybridization between two set
  of genotypes stored in \linkS4class{genind} objects (referred as the "2
  populations"). Allelic frequencies are derived for each population,
  and then gametes are sampled following a multinomial distribution. \cr

  The result consists in a set of 'n' genotypes, with different possible
  outputs (see 'res.type' argument).
}
\usage{
hybridize(x1, x2, n, res.type=c("genind","df","STRUCTURE"), file=NULL,
  quiet=FALSE, sep="/")
}
\arguments{
  \item{x1}{a \linkS4class{genind} object}
  \item{x2}{a \linkS4class{genind} object}
  \item{n}{an integer giving the number of hybrids requested}
  \item{res.type}{a character giving the type of output requested. Must
  be "genind" (default), "df" (i.e. data.frame like in
  \code{\link{genind2df}}), or "STRUCTURE" to generate a .str file
  readable by STRUCTURE (in which case the 'file' must be supplied). }
  \item{file}{a character giving the name of the file to be written
  when 'res.type' is "STRUCTURE"; if NULL, a the created file is of the
  form "hybrids\_[the current date].str".}
  \item{quiet}{a logical specifying whether the writing to a file (when
    'res.type' is "STRUCTURE") should be announced (FALSE, default) or
    not (TRUE).}
  \item{sep}{a character used to separate two alleles}
}
\value{
  A \linkS4class{genind} object (by default), or a data.frame of alleles
  (res.type="df"). No R output if res.type="STRUCTURE" (results written
  to the specified file).
  
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
## let's make a zebu-taurin hybrid
## P01 is Aubrac
## P15 is Zebu Peul
data(microbov)
microbov$pop.names

aubrac <- genind(truenames(microbov)$tab[microbov$pop=="P01",])
zebu <- genind(truenames(microbov)$tab[microbov$pop=="P15",])

## let's make some... Zebrac (or Aubbu)
zebrac <- hybridize(aubrac,zebu, n=40, res.type="df",sep="")
zebrac

## now let's merge all this in one data frame
newDat <- rbind.data.frame(genind2df(microbov), zebrac)
pop <- microbov$pop
levels(pop) <- microbov$pop.names
pop <- factor(c(as.character(pop),rep("zebrac",40)))

## make this dataset a genind to analyse it
obj <- df2genind(newDat, pop=pop)

## correspondance analysis
if(require(ade4)){
X <- genind2genpop(obj,missing="replace")
coa1 <- dudi.coa(as.data.frame(X$tab),scannf=FALSE,nf=3)
s.label(coa1$li,label=X$pop.names)
}
}
\keyword{manip}
