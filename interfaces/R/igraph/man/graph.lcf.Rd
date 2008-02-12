\name{graph.lcf}
\alias{graph.lcf}
\concept{LCF notation}
\title{Creating a graph from LCF notation}
\description{LCF is short for Lederberg-Coxeter-Frucht, it is a concise notation for
  3-regular Hamiltonian graphs. It constists of three parameters, the
  number of vertices in the graph, a list of shifts giving additional
  edges to a cycle backbone and another integer giving how many times
  the shifts should be performed. See
  \url{http://mathworld.wolfram.com/LCFNotation.html} for details.}
\usage{
graph.lcf(n, shifts, repeats)
}
\arguments{
  \item{n}{Integer, the number of vertices in the graph.}
  \item{shifts}{Integer vector, the shifts.}
  \item{repeats}{Integer constant, how many times to repeat the shifts.}
}
%\details{}
\value{A graph object.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{graph}} can create arbitrary graphs, see also the other
  functions on the its manual page for creating special graphs.
}
\examples{
# This is the Franklin graph:
g1 <- graph.lcf(12, c(5,-5), 6)
g2 <- graph.famous("Franklin")
graph.isomorphic.vf2(g1, g2)
}
\keyword{graphs}
