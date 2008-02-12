\name{graph.de.bruijn}
\alias{graph.de.bruijn}
\concept{De Bruijn graph}
\title{De Bruijn graphs.}
\description{De Bruijn graphs are labeled graphs representing the
  overlap of strings. }
\usage{
graph.de.bruijn(m,n)
}
\arguments{
  \item{m}{Integer scalar, the size of the alphabet. See details below.}
  \item{n}{Integer scalar, the length of the labels. See details below.}
}
\details{
 A de Bruijn graph represents relationships between strings. An alphabet
 of \code{m} letters are used and strings of length \code{n} are considered. 
 A vertex corresponds to every possible string and there is a directed edge
 from vertex \code{v} to vertex \code{w} if the string of \code{v} can
 be transformed into the string of \code{w} by removing its first letter
 and appending a letter to it.
 
 Please note that the graph will have \code{m} to the power \code{n}
 vertices and even more edges, so probably you don't want to supply too
 big numbers for \code{m} and \code{n}.
 
 De Bruijn graphs have some interesting properties, please see another source,
 eg. Wikipedia for details.   
}
\value{A graph object.}
\author{Gabor Csardi <csardi@rmki.kfki.hu>}
\seealso{\code{\link{graph.kautz}}, \code{\link{line.graph}}}
\examples{
# de Bruijn graphs can be created recursively by line graphs as well 
g <- graph.de.bruijn(2,1)
graph.de.bruijn(2,2)
line.graph(g)
}
\keyword{graphs}
