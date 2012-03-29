\name{graph.full.bipartite}
\alias{graph.full.bipartite}
\concept{Bipartite graph}
\concept{Two-mode network}
\title{Create a full bipartite graph}
\description{Bipartite graphs are also called two-mode by some. This
  function creates a bipartite graph in which every possible edge is
  present.}
\usage{
graph.full.bipartite (n1, n2, directed = FALSE, mode = c("all", "out", "in")) 
}
\arguments{
  \item{n1}{The number of vertices of the first kind.}
  \item{n2}{The number of vertices of the second kind.}
  \item{directed}{Logical scalar, whether the graphs is directed.}
  \item{mode}{Scalar giving the kind of edges to create for directed
    graphs. If this is \sQuote{\code{out}} then all vertices of the
    first kind are connected to the others; \sQuote{\code{in}} specifies
    the opposite direction; \sQuote{\code{all}} creates mutual
    edges. This argument is ignored for undirected graphs.x}
}
\details{
  Bipartite graphs have a \sQuote{\code{type}} vertex attribute in
  igraph, this is boolean and \code{FALSE} for the vertices of the first
  kind and \code{TRUE} for vertices of the second kind.
}
\value{
  An igraph graph, with the \sQuote{\code{type}} vertex attribute set.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.full}} for creating one-mode full graphs}
\examples{
g <- graph.full.bipartite(2, 3)
g2 <- graph.full.bipartite(2, 3, directed=TRUE)
g3 <- graph.full.bipartite(2, 3, directed=TRUE, mode="in")
g4 <- graph.full.bipartite(2, 3, directed=TRUE, mode="all")
}
\keyword{graphs}
