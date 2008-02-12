\name{get.adjlist}
\alias{get.adjlist}
\alias{get.adjedgelist}
\concept{Adjacency list}
\title{Adjacency lists}
\description{Create adjacency lists from a graph, either for adjacent
  edges or for neighboring vertices
}
\usage{
get.adjlist(graph, mode = c("all", "out", "in", "total"))
get.adjedgelist(graph, mode = c("all", "out", "in", "total"))
}
\arguments{
  \item{graph}{The input graph.}
  \item{mode}{Character scalar, it gives what kind of adjacent
    edges/vertices to include in the lists. \sQuote{\code{out}} is for
    outgoing edges/vertices, \sQuote{\code{in}} is for incoming
    edges/vertices, \sQuote{\code{all}} is for both. This argument is
    ignored for undirected graphs.}
}
\details{
  \code{get.adjlist} returns a list of numeric vectors, which include
  the ids of neighbor vertices (according to the \code{mode} argument) of all
  vertices.

  \code{get.adjedgelist} returns a list of numeric vectors, which
  include the ids of adjacent edgs (according to the \code{mode}
  argument) of all vertices.
}
\value{
  A list of numeric vectors.
}
%\references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{get.edgelist}}, \code{\link{get.adjacency}}}
\examples{
g <- graph.ring(10)
get.adjlist(g)
get.adjedgelist(g)
}
\keyword{graphs}
