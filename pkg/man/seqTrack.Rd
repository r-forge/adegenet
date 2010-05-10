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
  The SeqTrack algorithm aims at reconstructing genealogies of sampled
  haplotypes or genotypes for which a sampling date is
  available. Contrary to phylogenetics which aims at reconstructing
  hypothetical, unbserved most recent common ancestors (MRCA), SeqTrack
  considers that both ancestors and descendents are sampled together.

  This approach proved more efficient than phylogenetic approaches for
  reconstructing transmission trees in densely sampled disease outbreaks
  [1]. This implementation defines a generic function \code{seqTrack}

}
\usage{

\method{seqTrack}{matrix}(x, x.names, x.dates, best = c("min", "max"), prox.mat = NULL, 
    mu = NULL, haplo.length = NULL,  \dots)

plotSeqTrack(x, xy, useArrows = TRUE, annot = TRUE, labels = NULL, 
    dateRange = NULL, col = NULL, bg = "grey", add = FALSE, quiet = FALSE, 
    support = NULL, support.thres = 0.5, timeorder.thres = NULL, 
    mu = NULL, seq.length = NULL, col.pal = heat.colors, plot = TRUE, \dots)
}
\arguments{
  \item{seq.names}{}
  \item{seq.dates}{}
  \item{W}{}
  \item{optim}{}
  \item{x}{}
  \item{xy}{}
  \item{useArrows}{}
  \item{col}{}
  \item{bg}{}
  \item{add}{}
  \item{\dots}{}
  \item{}{}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\references{
  Jombart T, Eggo R, Dodd P, Balloux F (accepted) Reconstructing disease
  outbreaks from genetic data: a graph approach. Heredity.
}