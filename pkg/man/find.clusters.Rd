\encoding{UTF-8}
\name{find.clusters}
\alias{find.clusters}
\alias{find.clusters.data.frame}
\alias{find.clusters.matrix}
\alias{find.clusters.genind}
\title{}
\description{ == IN PROGRESS ==
  These functions implement the Discriminant Analysis of Principal Components
  (FIND.CLUSTERS). See 'details' section for a succint description of the method. FIND.CLUSTERS
  implementation calls upon \code{dudi.pca} from the \code{ade4} package and
  \code{lda} from the \code{MASS} package.

 \code{find.clusters} performs the FIND.CLUSTERS on a \code{data.frame}, a \code{matrix}, or a
 \code{\linkS4class{genind}} object, and returns an object with class
 \code{find.clusters}. If data are stored in a \code{data.frame} or a \code{matrix},
 these have to be quantitative data (i.e., \code{numeric} or \code{integers}),
 as opposed to \code{characters} or \code{factors}.

}
\usage{
\method{find.clusters}{data.frame}()

\method{find.clusters}{matrix}()

\method{find.clusters}{genind}()

}
\arguments{
\item{x}{\code{a data.frame}, \code{matrix}, or \code{\linkS4class{genind}}
  object. For the \code{data.frame} and \code{matrix} arguments, only
  quantitative variables should be provided.}
\item{grp,pop}{a \code{factor} indicating the group membership of individuals}
\item{n.pca}{an \code{integer} indicating the number of axes retained in the
  Principal Component Analysis (PCA) step. If \code{NULL}, interactive selection is triggered.}
\item{n.da}{an \code{integer} indicating the number of axes retained in the
  Discriminant Analysis step. If \code{NULL}, interactive selection is triggered.}
\item{center}{a \code{logical} indicating whether variables should be centred to
mean 0 (TRUE, default) or not (FALSE). Always TRUE for \linkS4class{genind} objects.}
\item{scale}{a \code{logical} indicating whether variables should be scaled
  (TRUE) or not (FALSE, default). Scaling consists in dividing variables by their
  (estimated) standard deviation to account for trivial differences in
  variances. Further scaling options are available for \linkS4class{genind}
  objects (see argument \code{scale.method}).}
\item{var.contrib,all.contrib}{a \code{logical} indicating whether the
  contribution of original variables (alleles, for \linkS4class{genind} objects)
  should be provided (TRUE) or not (FALSE, default). Such output can be useful,
  but can also create huge matrices when there the original size of the dataset
  is huge.}
\item{pca.select}{a \code{character} indicating the mode of selection of PCA
  axes, matching approximately "nbEig" or "percVar". For "nbEig", the user
  has to specify the number of axes retained (interactively, or via
  \code{n.pca}). For "percVar", the user has to specify the minimum amount of
  the total variance to be preserved by the retained axes, expressed as a
  percentage (interactively, or via \code{perc.pca}).  }
\item{perc.pca}{a \code{numeric} value between 0 and 100 indicating the
  minimal percentage of the total variance of the data to be expressed by the
  retained axes of PCA.}
\item{\ldots}{further arguments to be passed to other functions. For
  \code{find.clusters.matrix}, arguments are to match those of \code{find.clusters.data.frame}.}
\item{scale.method}{a \code{character} specifying the scaling method to be used
  for allele frequencies, which must match "sigma" (usual estimate of standard
  deviation) or "binom" (based on binomial distribution). See \code{\link{scaleGen}} for
  further details.}
\item{truenames}{a \code{logical} indicating whether true (i.e., user-specified)
  labels should be used in object outputs (TRUE, default) or not (FALSE).}
\item{xax,yax}{\code{integers} specifying which principal components of FIND.CLUSTERS
  should be shown in x and y axes. }
\item{col}{a suitable color to be used for groups. Not that the specified vector
should match the number of groups, not the number of individuals.}
\item{posi,bg,ratio,csub}{arguments used to customize the inset in scatterplots
  of FIND.CLUSTERS results. See \code{\link[pkg:ade4]{add.scatter}} documentation in the
  ade4 package for
  more details.}
\item{only.grp}{a \code{character} vector indicating which groups should be
  displayed. Values should match values of \code{x$grp}. If \code{NULL}, all
  results are displayed}
\item{subset}{\code{integer} or \code{logical} vector indicating which
  individuals should be displayed. If \code{NULL}, all
  results are displayed}
\item{cex.lab}{a \code{numeric} indicating the size of labels.}
\item{pch}{a \code{numeric} indicating the type of point to be used to indicate
  the prior group of individuals (see \code{\link{points}} documentation for
  more details).}
}
\details{

}
\value{
  The class \code{find.clusters} is a list with the following
  components:\cr
  \item{}{}

}
\references{
Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{PLoS genetics}.
}
\seealso{
    \code{\link{}}
    \code{\link{}}
    \code{\link{}}
    \code{\link{}}
    \code{\link{}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## data(find.clustersIllus), data(eHGDP), and data(H3N2) illustrate the find.clusters
## see ?find.clustersIllus, ?eHGDP, ?H3N2
##

example(find.clusters)


\dontrun{
example(eHGDP)
example(H3N2)
}

}
\keyword{multivariate}
