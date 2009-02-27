\name{page.rank}
\alias{page.rank}
\alias{page.rank.old}
\concept{Page rank}
\title{The Page Rank algorithm}
\description{ Calculates the Google PageRank for the specified
  vertices.}
\usage{
page.rank (graph, vids = V(graph), directed = TRUE, damping = 0.85, 
    weights = NULL, options = igraph.arpack.default) 
page.rank.old (graph, vids = V(graph), directed = TRUE, niter = 1000, 
    eps = 0.001, damping = 0.85, old = FALSE) 
}
\arguments{
  \item{graph}{The graph object. }
  \item{vids}{The vertices of interest.}
  \item{directed}{Logical, if true directed paths will be considered for
    directed graphs. It is ignored for undirected graphs.}
  \item{damping}{The damping factor (\sQuote{d} in the original paper).}
  \item{weights}{A numerical vector or \code{NULL}. This argument can be
    used to give edge weights for calculating the weighted PageRank of
    vertices. If this is \code{NULL} and the graph has a \code{weight}
    edge attribute then that is used. If \code{weights} is a numerical
    vector then it used, even if the graph has a \code{weights} edge
    attribute. If this is \code{NA}, then no edge weights are used (even
    if the graph has a \code{weight} edge attribute.}
  \item{options}{A named list, to override some ARPACK options. See
    \code{\link{arpack}} for details.}
  \item{niter}{The maximum number of iterations to perform.} 
  \item{eps}{The algorithm will consider the calculation as complete if
    the difference of PageRank values between iterations change less
    than this value for every node.}
  \item{old}{A logical scalar, whether the old style (pre igraph 0.5)
    normalization to use. See details below.}
}
\details{
  For the explanation of the PageRank algorithm, see the following
  webpage: \url{http://www-db.stanford.edu/~backrub/google.html}, or the
  following reference: 

  Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual
  Web Search Engine. Proceedings of the 7th World-Wide Web Conference,
  Brisbane, Australia, April 1998.  

  igraph 0.5 (and later) contains two PageRank calculation
  implementations. The \code{page.rank} function uses ARPACK to perform
  the calculation, see also \code{\link{arpack}}.

  The \code{page.rank.old} function performs a simple power method, this
  is the implementation that was available under the name
  \code{page.rank} in pre 0.5 igraph versions. Note that
  \code{page.rank.old} has an argument called \code{old}. If this
  argument is \code{FALSE} (the default), then the proper PageRank
  algorithm is used, i.e. \eqn{(1-d)/n} is added to the weighted
  PageRank of vertices to calculate the next iteration. If this
  argument is \code{TRUE} then \eqn{(1-d)} is added, just like in the
  PageRank paper; \eqn{d} is the damping factor, and \eqn{n} is the
  total number of vertices.
  A further difference is that the old implementation does not
  renormalize the page rank vector after each iteration.
  Note that the \code{old=FALSE} method is not stable,
  is does not necessarily converge to a fixed point. It should be
  avoided for new code, it is only included for compatibility with old
  igraph versions.
  
  Please note that the PageRank of a given vertex depends on the PageRank
  of all other vertices, so even if you want to calculate the PageRank for
  only some of the vertices, all of them must be calculated. Requesting
  the PageRank for only some of the vertices does not result in any
  performance increase at all. 

  Since the calculation is an iterative process, the algorithm is
  stopped after a given count of iterations or if the PageRank value
  differences between iterations are less than a predefined value. 
}
\value{
  For \code{page.rank} a named list with entries:
  \item{vector}{A numeric vector with the PageRank scores.}
  \item{value}{The eigenvalue corresponding to the eigenvector with
    the page rank scores. It should be always exactly one.}
  \item{options}{Some information about the underlying ARPACK
    calculation. See \code{\link{arpack}} for details.}
  
  For \code{page.rank.old} a numeric vector of Page Rank scores.
}
\references{
  Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual
  Web Search Engine. Proceedings of the 7th World-Wide Web Conference,
  Brisbane, Australia, April 1998.  
}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and
  Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{Other centrality scores: \code{\link{closeness}},
  \code{\link{betweenness}}, \code{\link{degree}}}
\examples{
g <- random.graph.game(20, 5/20, directed=TRUE)
page.rank(g)$vector

g2 <- graph.star(10)
page.rank(g2)$vector
}
\keyword{graphs}
