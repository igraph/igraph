\name{print.igraph}
\alias{print.igraph}
\alias{summary.igraph}
\title{Print graphs to the terminal}
\description{These functions attempt to print a graph to the terminal in
  a human readable form.}
\usage{
\method{print}{igraph}(x, \dots)
\method{summary}{igraph}(object, \dots)
}
\arguments{
  \item{x}{The graph to print.}
  \item{object}{The graph of which the summary will be printed.}
  \item{\dots}{Ignored.}
}
\details{
  \code{summary.igraph} prints the number of vertices, edges, whether the
  graph is directed and its storage type.
  \code{print.igraph} prints the same information as \code{summary.graph}
  and also the edges.
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
