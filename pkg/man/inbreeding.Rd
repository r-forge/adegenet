\encoding{UTF-8}
\name{Inbreeding}
\alias{inbreeding}
\title{Inbreeding coefficient for diploid genotypes}
\description{
  WARNING: this function is under development. Please contact the author
  (\email{t.jombart@imperial.ac.uk}) before using it.
  
  The function \code{inbreeding} computes Balloux's inbreeding
  coefficient for each individual of a \linkS4class{genind}
  objects. Results can be averaged over loci or detailed per locus. By
  default, \code{inbreeding} also produces a graphical output of the results.
}
\usage{
inbreeding(x, pop=NULL, truenames=TRUE, res.type=c("mean","byloc"), plot=TRUE, \ldots)
}
\arguments{
  \item{x}{an object of class \linkS4class{genind}.}
  \item{pop}{a factor giving the 'population' of each individual. If NULL,
    pop is seeked from \code{pop(x)}. Note that the term population refers in
    fact to any grouping of individuals'.}
  \item{truenames}{a logical indicating whether true names should be
    used (TRUE, default) instead of generic labels (FALSE); used if
    res.type is "matrix".}
 \item{res.type}{a character string matching "mean" or "byloc",
   specifying whether results should be averaged over loci ("mean") or
   detailed by locus ("byloc").}
 \item{plot}{a logical indicating whether a graphical
   output should be produced (TRUE, default), or not (FALSE).}
 \item{\ldots}{other arguments to be passed to \code{plot}.}
}
\value{
  A vector (if res.type is "mean"), or a matrix (if res.type is "byloc")
  of inbreeding coefficient values.
}
\seealso{
  \code{\link{Hs}}, \code{\link[hierfstat]{varcomp.glob}},
  \code{\link{gstat.randtest}}
}
\references{
  Brown AR, Hosken DJ, Balloux F, et al. 2009 Genetic variation,
  inbreeding and chemical exposure - combined effects in wildlife and
  critical considerations for ecotoxicology. Philosophical Transactions
  of the Royal Society B, London 364: 3377 - 3390
}
\details{
  Let \eqn{p_i} refer to the allele frequencies in a population. Let
  \eqn{h} be an variable which equates 1 if the individual is
  homozygote, and 0 otherwise. For one locus, Balloux's inbreeding coefficient is
  defined as:

  \eqn{  \frac{h - \sum_i p_i^2}{ \sum_i p_i^2 (1- \sum_i p_i^2)} } \cr

  For multi-locus genotypes, inbreeding values are averaged over the loci.
}
\author{
  Implementation: Thibaut Jombart \email{t.jombart@imperial.ac.uk}\cr
  Formula by Francois Balloux \email{f.balloux@imperial.ac.uk}
}
\examples{
## cat colonies of Nancy
data(nancycats)
inbreeding(nancycats)

## French/African cattle breeds
data(microbov)
inbreeding(microbov)
}
