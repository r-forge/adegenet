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
    \code{eHGDP} is a genind object with 3 supplementary components:
    \describe{
        \item{coun}{a factor giving the country of each individual (AF:
	  Afric; FR: France).}
        \item{breed}{a factor giving the breed of each individual.}
        \item{spe}{is a factor giving the species of each individual
	  (BT: Bos taurus; BI: Bos indicus).}
    }
}
\source{
Data prepared by Katayoun Moazami-Goudarzi and Denis Lalo\"e (INRA,
Jouy-en-Josas, France)
}
\references{
  Lalo\"e D., Jombart T., Dufour A.-B. and Moazami-Goudarzi K. (2007)
  Consensus genetic structuring and typological value of markers using
  Multiple Co-Inertia Analysis. \emph{Genetics Selection Evolution}.
  \bold{39}: 545--567.
}
\examples{
data(eHGDP)
eHGDP



}
\keyword{datasets}
