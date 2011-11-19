\name{scgExtra}
\alias{scg.normEps}
\title{SCG Extra Functions}
\description{
  Some useful functions to perform general actions in Spectral Coarse
  Graining (SCG).
}
\usage{
scg.normEps(V, groups, matrix = c("symmetric", "laplacian",
          "stochastic"), p = NULL, norm = c("row", "col"))
}
\arguments{
  \item{V}{A numeric matrix of (eigen)vectors assumed normalized.
    The vectors are to be stored column-wise in \code{V}).}
  \item{groups}{A vector of \code{nrow(V)} integers labeling 
    each group vertex in the partition.}
  \item{matrix}{The type of semi-projector used for the SCG. For 
    now \dQuote{symmetric}, \dQuote{laplacian} and \dQuote{stochastic}
    are available.} 
  \item{p}{A probability vector of length \code{nrow(V)}.
    \code{p} is the stationary probability distribution of a Markov chain 
    when \code{matrix} = \dQuote{stochastic}. This parameter is ignored otherwise.}
  \item{norm}{Either \dQuote{row} or \dQuote{col}. If set to \dQuote{row} 
    the rows of the Laplacian matrix sum to zero and the rows of the 
    stochastic matrix sum to one; otherwise it is the columns.}
}
\details{
 \code{scg.normEps} computes \eqn{\Vert v_i-Pv_i\Vert}, where \eqn{v_i}
  is the \eqn{i}th eigenvector in \code{V} and \eqn{P} is the projector
  corresponding to the \code{matrix} argument.
}
\value{
  \code{normEps} returns with a numeric vector whose \eqn{i}th component is 
  \eqn{\Vert v_i-Pv_i\Vert} (see Details).
}
\references{
  D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios, Shrinking
  Matrices while Preserving their Eigenpairs with Application to the
  Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on
    Matrix Analysis and Applications}, 2008. 
  \url{http://people.epfl.ch/david.morton}
}
\author{David Morton de Lachapelle \email{david.morton@epfl.ch},
  \email{david.mortondelachapelle@swissquote.ch}}
\seealso{\link{SCG} and \code{\link{scg}}. }
\examples{
v <- rexp(20)
km <- kmeans(v,5)
sum(km$withinss)
normEps(v,km$cluster)^2
}
\keyword{graphs}
