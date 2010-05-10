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
\alias{get.likelihood.seqTrack.default}
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

plotSeqTrack(x, xy, useArrows = TRUE, annot = TRUE, labels = NULL, 
    dateRange = NULL, col = NULL, bg = "grey", add = FALSE, quiet = FALSE, 
    support = NULL, support.thres = 0.5, timeorder.thres = NULL, 
    mu = NULL, seq.length = NULL, col.pal = heat.colors, plot = TRUE, \dots)
}
\arguments{
  \item{x}{for seqTrack, a matrix giving weights to pairs of ancestries such that
    x[i,j] is the weight of 'i ancestor of j'. For plotSeqTrack, a
    seqTrack object.}
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
  \item{useArrows}{a logical indicating whether arrows should be used to
    represented ancestries (pointing from ancestor to descendent, TRUE),
    or whether segments shall be used (FALSE).}
  \item{annot}{a logical indicating whether arrows or segments
    representing ancestries should be annotated (TRUE) or not (FALSE).}
  \item{labels}{a character vector containing annotations of the
    ancestries. If left empty, ancestries are annotated by the descendent.}
  \item{dateRange}{a vector of length two with POSIXct format indicating
    the time window for which ancestries should be displayed.}
  \item{col}{a vector of colors to be used for plotting ancestries.}
  \item{bg}{a color to be used as background.}
  \item{add}{a logical stating whether the plot should be added to
    current figure (TRUE), or drawn as a new plot (FALSE, default).}
  \item{quiet}{a logical stating whether messages other than errors
    should be displayed (FALSE, default), or hidden (TRUE).}
  \item{support}{an optional vector of numbers between 0 and 1
    indicating the support for each ancestry. If provided and if
    \code{col} is left empty, ancestries are represented with colors
    indicating the support.}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  
  \item{col}{}
  \item{bg}{}
  \item{add}{}
  \item{\dots}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\references{
  Jombart T, Eggo R, Dodd P, Balloux F (accepted) Reconstructing disease
  outbreaks from genetic data: a graph approach. Heredity.
}