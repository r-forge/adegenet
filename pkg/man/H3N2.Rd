\encoding{UTF-8}
\name{H3N2}
\alias{H3N2}
\docType{data}
\title{}
\description{
This dataset consists of 2949 strains of seasonal influenza (H3N2)
distributed worldwide, and typed at 125 SNPs located in the
hemagglutinin (HA) segment. These data were gathered from DNA sequences
available from Genbank (http://www.ncbi.nlm.nih.gov/Genbank/).
}
\usage{data(H3N2)}
\format{
    \code{H3N2} is a genind object with a data frame named \code{...} as supplementary
    component (\code{H3N2@other)}, which contains the following items:
    \describe{
        \item{x}{a \code{data.frame} containing miscellanous annotations
	  of the sequences.}
        \item{Region}{a character vector indicating the geographic region of each population.}
        \item{Label}{a character vector indicating the correspondance
	  with population labels used in the genind object (i.e., as
	  output by \code{pop(H3N2)}).}
	\item{Latitude,Longitude}{geographic coordinates of the
	  populations, indicated as north and east degrees.}
    }
}
\source{
  This dataset was prepared by Thibaut Jombart
  (t.jombart@imperia.ac.uk), from annotated sequences available on
  Genbank (http://www.ncbi.nlm.nih.gov/Genbank/).
}
\references{
Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{PLoS genetics}.
}
\examples{
\dontrun{
## LOAD DATA
data(H3N2)
H3N2


## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc1 <- dapc(H3N2, all.contrib=TRUE, scale=FALSE, n.pca=200, n.da=80) # takes 2 minutes
dapc1

## (see ?dapc for details about the output)



## SCREEPLOT OF EIGENVALUES
barplot(dapc1$eig, main="H3N2 - DAPC eigenvalues", col=c("red","green","blue", rep("grey", 1000)))



## SCATTERPLOTS
## (!) Note: colors may be inverted with respect to the
## original paper (as signs of principal components are arbitrary)
## axes 1-2
s.label(dapc1$grp.coord[,1:2], clab=0, sub="Axes 1-2")
par(xpd=T)
colorplot(dapc1$grp.coord[,1:2], dapc1$grp.coord, cex=3, add=TRUE)
add.scatter.eig(dapc1$eig,10,1,2, posi="bottomright", ratio=.3, csub=1.25)

## axes 2-3
s.label(dapc1$grp.coord[,2:3], clab=0, sub="Axes 2-3")
par(xpd=T)
colorplot(dapc1$grp.coord[,2:3], dapc1$grp.coord, cex=3, add=TRUE)
add.scatter.eig(dapc1$eig,10,1,2, posi="bottomright", ratio=.3, csub=1.25)



## MAP DAPC1 RESULTS
if(require(maps)){

xy <- cbind(H3N2$other$popInfo$Longitude, H3N2$other$popInfo$Latitude)

par(mar=rep(.1,4))
map(fill=TRUE, col="lightgrey")
colorplot(xy, -dapc1$grp.coord, cex=3, add=TRUE, trans=FALSE)
}



## LOOK FOR OTHER CLUSTERS
## to reproduce results of the reference paper, use :
## grp <- find.clusters(hgdp, max.n=50, n.pca=200, scale=FALSE)
## and then
## plot(grp$Kstat, type="b", col="blue")

grp <- find.clusters(H3N2, max.n=30, n.pca=200, scale=FALSE, n.clust=4) # takes about 2 minutes
names(grp)

## (see ?find.clusters for details about the output)



## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc2 <- dapc(H3N2, pop=grp$grp, all.contrib=TRUE, scale=FALSE, n.pca=200, n.da=80) # takes around 2 minutes
dapc2


## PRODUCE SCATTERPLOT
scatter(dapc2) # axes 1-2
scatter(dapc2,2,3) # axes 2-3


## MAP DAPC2 RESULTS
if(require(maps)){
xy <- cbind(H3N2$other$popInfo$Longitude, H3N2$other$popInfo$Latitude)

myCoords <- apply(dapc2$ind.coord, 2, tapply, pop(H3N2), mean)

par(mar=rep(.1,4))
map(fill=TRUE, col="lightgrey")
colorplot(xy, myCoords, cex=3, add=TRUE, trans=FALSE)
}

}
}
\keyword{datasets}
