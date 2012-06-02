\name{print.igraph}
\alias{print.igraph}
\alias{str.igraph}
\alias{summary.igraph}
\title{Print graphs to the terminal}
\description{These functions attempt to print a graph to the terminal in
  a human readable form.}
\usage{
\method{print}{igraph}(x, full=getIgraphOpt("print.full"),
  graph.attributes=getIgraphOpt("print.graph.attributes"),
  vertex.attributes=getIgraphOpt("print.vertex.attributes"),
  edge.attributes=getIgraphOpt("print.edge.attributes"),
  names=TRUE, \dots)
\method{summary}{igraph}(object, \dots)
\method{str}{igraph}(object, \dots)
}
\arguments{
  \item{x}{The graph to print.}
  \item{full}{Logical scalar, whether to print the graph structure
    itself as well.}
  \item{graph.attributes}{Logical constant, whether to print graph attributes.}
  \item{vertex.attributes}{Logical constant, whether to print vertex attributes.}
  \item{edge.attributes}{Logical constant, whether to print edge
    attributes.}
  \item{names}{Logical constant, whether to print symbolic vertex names
    (ie. the \code{name} vertex attribute) or vertex ids.}
  \item{object}{The graph of which the summary will be printed.}
  \item{\dots}{Additional agruments.}
}
\details{
  \code{summary.igraph} prints the number of vertices, edges and whether the
  graph is directed.
  \code{str.igraph} prints the same information, and also lists
  the edges, and optionally graph, vertex and/or edge attributes.
  \code{print.igraph} behaves either as \code{summary.igraph} or
  \code{str.igraph} depending on the \code{full} argument. See also the
  \sQuote{print.full} igraph option and \code{\link{getIgraphOpt}}.
  
  As of igraph 0.4 \code{print.igraph} uses the \code{max.print} option,
  see \code{\link[base]{options}} for details.
}
\value{All these functions return the graph invisibly.}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
%\seealso{}
\examples{
g <- graph.ring(10)
g
summary(g)
}
\keyword{graphs}
