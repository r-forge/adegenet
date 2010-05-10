\name{seqTrack}
\alias{seqTrack}
% \alias{seqTrackG}
\alias{seqTrack.default}
% \alias{seqTrackG.default}
\alias{plotSeqTrack}
% \alias{optimize.seqTrack}
% \alias{optimize.seqTrack.default}
\alias{get.likelihood}
\alias{get.likelihood.seqTrack}
\alias{seqTrack-class}
% \alias{get.result.by}
% \alias{get.consensus}
\alias{as,seqTrack,graphNEL-method}
\alias{coerce,seqTrack,graphNEL-method}
\title{SeqTrack algorithm for reconstructing genealogies}
\description{
  The SeqTrack algorithm [1] aims at reconstructing genealogies of sampled
  haplotypes or genotypes for which a collection date is
  available. Contrary to phylogenetic methods which aims at reconstructing
  hypothetical ancestors for observed sequences, SeqTrack
  considers that ancestors and descendents are sampled together, and
  therefore infers ancestry relationships among the sampled sequences.

  This approach proved more efficient than phylogenetic approaches for
  reconstructing transmission trees in densely sampled disease outbreaks
  [1]. This implementation defines a generic function \code{seqTrack}
  with methods for specific object classes.

}
\usage{
seqTrack(\ldots)

\method{seqTrack}{matrix}(x, x.names, x.dates, best = c("min", "max"), prox.mat = NULL, 
    mu = NULL, haplo.length = NULL,  \dots)

plotSeqTrack(x, xy, use.arrows=TRUE, annot=TRUE, labels=NULL,
                         col=NULL, bg="grey", add=FALSE, quiet=FALSE, date.range=NULL,
                         plot=TRUE, \dots)
}
\arguments{
  \item{x}{for seqTrack, a matrix giving weights to pairs of ancestries
    such that x[i,j] is the weight of 'i ancestor of j'. For
    plotSeqTrack, a seqTrack object.}
  \item{seq.names}{a character vector giving the labels of the haplotypes/genotypes}
  \item{seq.dates}{a vector of collection dates for the sampled
    haplotypes/genotypes. Dates must have the POSIXct format. See
    \code{details} or \code{?as.POSIXct} for more information.}
  \item{best}{a character string matching 'min' or 'max', indicating
    whether genealogies should minimize or maximize the sum of weights
    of ancestries.}
  \item{prox.mat}{an optional matrix of proximities between
    haplotypes/genotypes used to resolve ties in the choice of
    ancestors, by picking up the 'closest' ancestor amongst possible
    ancestors, in the sense of \code{prox.mat}. \code{prox.mat[i,j]}
    must indicate a proximity for the relationship 'i ancestor to
    j'. For instance, if \code{prox.mat} contains spatial proximities,
    then \code{prox.mat[i,j]} gives a measure of how easy it is to
    migrate from location 'i' to 'j'.
  }
  \item{mu}{(optional) a mutation rate, per site and per day. When 'x'
    contains numbers of mutations, used to resolve ties using a maximum
    likelihood approach (requires \code{haplo.length} to be provided).}
 \item{haplo.length}{(optional) the length of analysed sequences in
    number of nucleotides. When 'x' contains numbers of mutations, used
    to resolve ties using a maximum likelihood approach (requires
    \code{mu} to be provided).}  
  \item{xy}{spatial coordinates of the sampled haplotypes/genotypes.}
  \item{use.arrows}{a logical indicating whether arrows should be used to
    represented ancestries (pointing from ancestor to descendent, TRUE),
    or whether segments shall be used (FALSE).}
  \item{annot}{a logical indicating whether arrows or segments
    representing ancestries should be annotated (TRUE) or not (FALSE).}
  \item{labels}{a character vector containing annotations of the
    ancestries. If left empty, ancestries are annotated by the descendent.}
  \item{col}{a vector of colors to be used for plotting ancestries.}
  \item{bg}{a color to be used as background.}
  \item{add}{a logical stating whether the plot should be added to
    current figure (TRUE), or drawn as a new plot (FALSE, default).}
  \item{quiet}{a logical stating whether messages other than errors
    should be displayed (FALSE, default), or hidden (TRUE).}
   \item{date.range}{a vector of length two with POSIXct format indicating
    the time window for which ancestries should be displayed.}
  \item{plot}{a logical stating whether a plot should be drawn (TRUE,
    default), or not (FALSE). In all cases, the
    function invisibly returns plotting information.}
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\references{
  Jombart T, Eggo R, Dodd P, Balloux F (accepted) Reconstructing disease
  outbreaks from genetic data: a graph approach. Heredity.
}
\value{
  === output of seqTrack ===
  seqTrack function returns data.frame with the class \code{seqTrack},
  in which each row is an inferred ancestry described by the following columns:
  - id: indices identifying haplotypes/genotypes
  - ances: index of the inferred ancestor
  - weight: weight of the inferred ancestries
  - date: date of the haplotype/genotype
  - ances.date: date of the ancestor

  === output of plotSeqTrack ===
  This graphical function invisibly returns the coordinates of the
  arrows/segments drawn and their colors, as a data.frame.
}
\details{
  === Maximum parsimony genealogies ===
  
  === Likelihood of observed genetic differentiation ===
  === Converting seqTrack objects to graphs ===

  
}
\seealso{
  \code{\link{get.likelihood.seqTrack}}
}