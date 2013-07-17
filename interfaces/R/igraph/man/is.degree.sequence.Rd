\name{is.degree.sequence}
\alias{is.degree.sequence}
\alias{is.graphical.degree.sequence}
\concept{Degree sequence}
\title{Degree sequences of graphs}
\description{
  These functions decide whether a sequence (or two, if the graph is
  directed) of integers can be realized as vertex degrees by a graph or
  simple graph.
}
\usage{
is.degree.sequence (out.deg, in.deg = NULL)
is.graphical.degree.sequence (out.deg, in.deg = NULL)
}
\arguments{
  \item{out.deg}{Integer vector, the degree sequence for undirected
    graphs, or the out-degree sequence for directed graphs.}
  \item{in.deg}{\code{NULL} or an integer vector. For undireted graphs,
    it should be \code{NULL}. For directed graphs it specifies the
    in-degrees.}
}
\details{
  \code{is.degree.sequence} checks whether the given vertex degrees (in-
  and out-degrees for directed graphs) can be realized by a graph. Note
  that the graph does not have to be simple, it may contain loop and
  multiple edges. For undirected graphs, it also checks whether the sum
  of degrees is even. For directed graphs, the function checks whether
  the lengths of the two degree vectors are equal and whether their sums
  are also equal. These are known sufficient and necessary conditions
  for a degree sequence to be valid.

  \code{is.graphial.degree.sequence} determines whether the given vertex
  degrees (in- and out-degrees for directed graphs) can be reliazed in a
  simple graph, i.e. a graph without multiple or loop edges.
}
\value{
  A logical scalar.
}
\references{
  Hakimi SL: On the realizability of a set of integers as degrees of the
  vertices of a simple graph. \emph{J SIAM Appl Math} 10:496-506, 1962.
 
  PL Erdos, I Miklos and Z Toroczkai: A simple Havel-Hakimi type algorithm
  to realize graphical degree sequences of directed graphs.
  \emph{The Electronic Journal of Combinatorics} 17(1):R66, 2010.
}
\author{
  Tamas Nepusz \email{ntamas@gmail.com}
}
\examples{
g <- erdos.renyi.game(100, 2/100)
is.degree.sequence(degree(g))
is.graphical.degree.sequence(degree(g))
}
\keyword{graphs}
