\name{print.graph}
\alias{print.graph}
\alias{summary.graph}
\title{Print graphs to the terminal}
\description{These functions attempt to print a graph to the terminal in
  a human readable form.}
\usage{
\method{print}{graph}(graph)
\method{summary}{graph}(graph)
}
\arguments{
  \item{graph}{The graph to print.}
}
\details{
  \code{summary.graph} prints the number of vertices, edges, whether the
  graph is directed and its storage type.
  \code{print.graph} prints the same information as \code{summary.graph}
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
