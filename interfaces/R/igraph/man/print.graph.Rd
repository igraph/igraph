\name{print.igraph}
\alias{print.igraph}
\alias{summary.igraph}
\title{Print graphs to the terminal}
\description{These functions attempt to print a graph to the terminal in
  a human readable form.}
\usage{
\method{print}{igraph}(x, graph.attributes=FALSE,
  vertex.attributes=FALSE, edge.attributes=FALSE, \dots)
\method{summary}{igraph}(object, \dots)
}
\arguments{
  \item{x}{The graph to print.}
  \item{graph.attributes}{Logical constant, whether to print graph attributes.}
  \item{vertex.attributes}{Logical constant, whether to print vertex attributes.}
  \item{edge.attributes}{Logical constant, whether to print edge attributes.}
  \item{object}{The graph of which the summary will be printed.}
  \item{\dots}{Ignored.}
}
\details{
  \code{summary.igraph} prints the number of vertices, edges and whether the
  graph is directed.
  \code{print.igraph} prints the same information, aand also lists
  the edges, and optionally graph, vertex and/or edge attributes.
}
\value{Both functions return the graph invisibly.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
%\seealso{}
\examples{
g <- graph.ring(10)
g
summary(g)
}
\keyword{graphs}
