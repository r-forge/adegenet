\encoding{UTF-8}
\name{eHGDP}
\alias{eHGDP}
\docType{data}
\title{Microsatellites genotypes of 15 cattle breeds}
\description{
This data set gives the genotypes of 704 cattle individuals for 30
microsatellites recommended by the FAO. The individuals are divided into
two countries (Afric, France), two species (Bos taurus, Bos indicus) and
15 breeds. Individuals were chosen in order to avoid pseudoreplication
according to their exact genealogy.
}
\usage{data(eHGDP)}
\format{
    \code{eHGDP} is a genind object with a data frame named \code{popInfo} as supplementary
    component (\code{eHGDP@other$popInfo)}, which contains the following
    variables:
    \describe{
        \item{Population}{a character vector indicating populations.}
        \item{Region}{a character vector indicating the geographic region of each population.}
        \item{Label}{a character vector indicating the correspondance
	  with population labels used in the genind object (i.e., as
	  output by \code{pop(eHGDP)}).}
	\item{Latitude,Longitude}{geographic coordinates of the
	  populations, indicated as north and east degrees.}
    }
}
\source{
Data prepared by Francois Balloux.
}
\references{
  Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{PLoS genetics}.
}
\examples{
\dontrun{
## LOAD DATA
data(eHGDP)
eHGDP


## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc1 <- dapc(eHGDP, all.contrib=TRUE, scale=FALSE, n.pca=200, n.da=80) # takes 2 minutes
dapc1

## SCREEPLOT OF EIGENVALUES
barplot(dapc1$eig, main="eHGDP - DAPC eigenvalues",
col=c("red","green","blue", rep("grey", 1000)))

## SCATTERPLOTS
## ! note ! colors may be inverted with respect to the
## original paper (signs of principal components are arbitrary)
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


## MAP DAPC RESULTS
if(require(maps)){

xy <- cbind(eHGDP$other$popInfo$Longitude, eHGDP$other$popInfo$Latitude)

par(mar=rep(.1,4))
map(fill=TRUE, col="lightgrey")
colorplot(xy, dapc1$grp.coord, cex=3, add=TRUE, trans=FALSE)
}



## LOOK FOR LARGER CLUSTERS
## to reproduce results of the reference paper, use :
## grp <- find.clusters(hgdp, max.n=50, n.pca=200, scale=FALSE)
## and then
## plot(grp$Kstat, type="b", col="blue")

grp <- find.clusters(hgdp, max.n=30, n.pca=200, scale=FALSE, n.clust=4) # takes about 2 minutes
names(grp)


## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc1 <- dapc(eHGDP, grp=grp, all.contrib=TRUE, scale=FALSE, n.pca=200, n.da=80) # takes around 2 minutes
dapc1



}
}
\keyword{datasets}
