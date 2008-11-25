\name{scgExtra}
\alias{normEps}
\alias{groupComposition}
\alias{renumber}
\title{SCG Extra Functions}
\description{
Some useful functions to perform general actions in Spectral Coarse Graining (SCG). 
}

\usage{
normEps(V,  gr,  matrix = c("symmetric", "laplacian", "stochastic"),
	norm = "row",  p = NULL)
groupComposition(gr, vertex.smallest = 1)
renumber(x, smallest = 1, order = FALSE)
}

\arguments{
  \item{V}{A numeric matrix of (eigen)vectors assumed normalized.
  The vectors are to be stored column-wise in \code{V}).}

  \item{gr}{A vector of \code{nrow(V)} integers labeling 
  each group vertex in the partition.}

  \item{matrix}{The type of semi-projector used for the SCG. For 
  now \dQuote{symmetric}, \dQuote{laplacian} and \dQuote{stochastic} are available.}

  \item{norm}{Either \dQuote{row} or \dQuote{col}. If set to \dQuote{row} 
  the rows of the Laplacian matrix sum to zero and the rows of the 
  stochastic matrix sum to one; otherwise it is the columns.}

  \item{p}{A probability vector of length \code{nrow(V)}.
  \code{p} is the stationary probability distribution of a Markov chain 
  when \code{matrix} = \dQuote{stochastic}. This parameter is ignored otherwise.}

  \item{vertex.smallest}{An integer used for the smallest vertex index.}

  \item{x}{A vector of labels allowed to be cast as \code{factor}.}

  \item{smallest}{The smallest integer in the renumbered vector.}

  \item{order}{If \code{TRUE} the renumbered vector will be sorted in increasing order 
  (with \code{smallest} in the first position).}
}

\details{
 \code{normEps} computes \eqn{\Vert v_i-Pv_i\Vert}, where \eqn{v_i} is the
 \eqn{i}th eigenvector in \code{V} and \eqn{P} is the projector
  corresponding to the \code{matrix} argument.
  
  \code{renumber} Applies successively \code{as.factor} and \code{as.integer} 
  to \code{x}. Hence, labels of any type that can be cast as \code{factor} are allowed 
  in \code{x} (e.g. \code{integer} or \code{character}).
}
\value{
  \code{normEps} returns with a numeric vector whose \eqn{i}th component is 
  \eqn{\Vert v_i-Pv_i\Vert} (see Details).
  \code{groupComposition} returns with the list of vertices in each group.
  \code{renumber} returns with a vector of positive integers.
}

\references{
D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios, Shrinking Matrices while Preserving their Eigenpairs with Application to the Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on Matrix Analysis and Applications}, 2008.
\url{http://people.epfl.ch/david.morton}
}
\author{David Morton de Lachapelle \email{david.morton@epfl.ch}, \email{david.mortondelachapelle@swissquote.ch}}
\seealso{\code{\link{scg}} }

\examples{
v <- rexp(20)
km <- kmeans(v,5)
sum(km$withinss)
normEps(v,km$cluster)^2

km$cluster
groupComposition(km$cluster)

renumber(c("hi",3,"hi",-4,1,3), order=TRUE)
}

\keyword{array}
\keyword{graphs}
