\name{page.rank}
\alias{page.rank}
\title{The Page Rank algorithm}
\description{ Calculates the Google PageRank for the specified
  vertices.}
\usage{
page.rank(graph, nodes = V(graph), directed = is.directed(graph),
    niter = 1000, eps = 0.001, damping = 0.85)
}
\arguments{
  \item{graph}{The graph object. }
  \item{nodes}{The vertices of interest.}
  \item{directed}{Logical, if true directed paths will be considered for
    directed graphs. It is ignored for undirected graphs.}
  \item{niter}{The maximum number of iterations to perform.} 
  \item{eps}{The algorithm will consider the calculation as complete if
    the difference of PageRank values between iterations change less
    than this value for every node.}
  \item{damping}{The damping factor (\sQuote{d} in the original paper).}
}
\details{
  Please note that the PageRank of a given vertex depends on the PageRank
  of all other vertices, so even if you want to calculate the PageRank for
  only some of the vertices, all of them must be calculated. Requesting
  the PageRank for only some of the vertices does not result in any
  performance increase at all. 

  Since the calculation is an iterative process, the algorithm is
  stopped after a given count of iterations or if the PageRank value
  differences between iterations are less than a predefined value. 

  For the explanation of the PageRank algorithm, see the following
  webpage: \url{http://www-db.stanford.edu/~backrub/google.html}, or the
  following reference: 

  Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual
  Web Search Engine. Proceedings of the 7th World-Wide Web Conference,
  Brisbane, Australia, April 1998.  
}
\value{A numeric vector of Page Rank scores.}
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
page.rank(g)

g2 <- graph.star(10)
page.rank(g2)
}
\keyword{graphs}

