\name{graph.kautz}
\alias{graph.kautz}
\concept{Kautz graph}
\title{Kautz graphs}
\description{Kautz graphs are labeled graphs representing the
  overlap of strings. }
\usage{
graph.kautz(m,n)
}
\arguments{
  \item{m}{Integer scalar, the size of the alphabet. See details below.}
  \item{n}{Integer scalar, the length of the labels. See details below.}
}
\details{
  A Kautz graph is a labeled graph, vertices are labeled by strings
  of length \code{n+1} above an alphabet with \code{m+1} letters, with
  the restriction that every two consecutive letters in the string
  must be different. There is a directed edge from a vertex \code{v} to
  another vertex \code{w} if it is possible to transform the string of
  \code{v} into the string of \code{w} by removing the first letter and
  appending a letter to it.
  
  Kautz graphs have some interesting properties, see eg. Wikipedia
  for details.
}
\value{A graph object.}
\author{Gabor Csardi <csardi@rmki.kfki.hu>, the first version in R was
  written by Vincent Matossian.}
\seealso{\code{\link{graph.de.bruijn}}, \code{\link{line.graph}}}
\examples{
line.graph(graph.kautz(2,1))
graph.kautz(2,2)
}
\keyword{graphs}
