\encoding{UTF-8}
\name{findMutations}
\alias{findMutations}
\alias{findMutations.DNAbin}
\title{Identify mutations between DNA sequences}
\description{
  This function compares pairs of aligned DNA sequences and identify
  mutations (position and nature).\cr

  The function \code{findMutations} is a generic, but the only method
  implemented in adegenet so far is for \code{\link[ape]{DNAbin}}
  objects.
}
\usage{
findMutations(x, pairs=NULL, \dots)

\method{findMutations}{DNAbin}(x, pairs=NULL, \dots)

}
\arguments{
  \item{x}{a \code{DNAbin} object containing aligned sequences, as a matrix.}
  \item{pairs}{a matrix with two columns specifying which
    sequences to compare. Mutations will be identified from the first
    column to the second. Any way of refering to the rows of \code{x} is
    acceptable. If \code{NULL}, all pairs are considered.}
}
\value{
  A named list indicating the mutations from one sequence to another.
  For each comparison, a three-column matrix is provided, corresponding
  to the nucleotides in first and second sequence, and a summary of the
  mutation provided as:
  [position]:[nucleotide in first sequence]->[nucleotide in second sequence].
}
\seealso{
  The \code{\link{fasta2DNAbin}} to read fasta alignments with minimum
  RAM use.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}.
 }
\examples{
\dontrun{
if(require(ape)){
data(woodmouse)
findMutations(woodmouse[1:5,])
findMutations(woodmouse, pairs=cbind(c(1,1),c(2,3)))
}
}
}
