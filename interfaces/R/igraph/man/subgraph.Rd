\name{subgraph}
\alias{subgraph}
\alias{induced.subgraph}
\concept{Subgraph}
\title{Subgraph of a graph}
\description{\code{subgraph} creates a subgraph of a graph, containing
  only the specified vertices and all the edges among them.}
\usage{
induced.subgraph(graph, vids, impl=c("auto", "copy_and_delete",
     "create_from_scratch")
subgraph(graph, v)
}
\arguments{
  \item{graph}{The original graph.}
  \item{vids,v}{Numeric vector, the vertices of the original 
    graph which will form the subgraph.}
  \item{impl}{Character scalar, to choose between two 
    implementation of the subgraph
    calculation. \sQuote{\code{copy_and_delete}} copies the graph
    first, and then deletes the vertices and edges that are not
    included in the result graph. \sQuote{\code{create_from_scratch}}
    searches for all vertices and edges that must be kept and then
    uses them to create the graph from scratch. \sQuote{\code{auto}}
    chooses between the two implementations automatically, using
    heuristics based on the size of the original and the result
    graph.}
}
\details{
  \code{induced.subgraph} calculates the induced subgraph of a set of
  vertices in a graph. This means that exactly the specified vertices
  and all the edges between then will be kept in the result graph.

  The \code{subgraph} function does the same currently, but it is
  deprecated and will be removed in the next major version of igraph.
}
\value{A new graph object.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso
\examples{
g <- graph.ring(10)
g2 <- induced.subgraph(g, 1:7)
}
\keyword{graphs}
