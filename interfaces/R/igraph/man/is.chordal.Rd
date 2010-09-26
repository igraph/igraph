\name{is.chordal}
\alias{is.chordal}
\concept{maximum cardinality search}
\concept{graph decomposition}
\concept{chordal graph}
\title{Chordality of a graph}
\description{
  A graph is chordal (or triangulated) if each of its cycles of four or
  more nodes has a chord, which is an edge joining two nodes that are
  not adjacent in the cycle. An equivalent definition is that any
  chordless cycles have at most three nodes.}
\usage{
is.chordal(graph, alpha = NULL, alpham1 = NULL, fillin = FALSE,
           newgraph = FALSE)
}
\arguments{
  \item{graph}{The input graph. It may be directed, but edge directions
    are ignored, as the algorithm is defined for undirected graphs.}
  \item{alpha}{Numeric vector, the maximal chardinality ordering of the
    vertices. If it is \code{NULL}, then it is automatically calculated
    by calling \code{\link{maximum.cardinality.search}}, or from
    \code{alpham1} if that is given..}
  \item{alpham1}{Numeric vector, the inverse of \code{alpha}. If it is
    \code{NULL}, then it is automatically calculated by calling
    \code{\link{maximum.cardinality.search}}, or from \code{alpha}.}
  \item{fillin}{Logical scalar, whether to calculate the fill-in edges.}
  \item{newgraph}{Logical scalar, whether to calculate the triangulated
    graph.}
}
\details{
  The chordality of the graph is decided by first performing maximum
  cardinality search on it (if the \code{alpha} and \code{alpham1}
  arguments are \code{NULL}), and then calculating the set of fill-in
  edges.

  The set of fill-in edges is empty if and only if the graph is
  chordal.

  It is also true that adding the fill-in edges to the graph makes it
  chordal.
}
\value{
  A list with three members:
  \item{chordal}{Logical scalar, it is \code{TRUE} iff the input graph
    is chordal.}
  \item{fillin}{If requested, then a numeric vector giving the fill-in
    edges. \code{NULL} otherwise.}
  \item{newgraph}{If requested, then the triangulated graph, an
    \code{igraph} object. \code{NULL} otherwise.}
}
\references{ Robert E Tarjan and Mihalis Yannakakis. (1984). Simple
  linear-time algorithms to test chordality of graphs, test acyclicity
  of hypergraphs, and selectively reduce acyclic hypergraphs.
  \emph{SIAM Journal of Computation} 13, 566--579.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{maximum.cardinality.search}} }
\examples{
## The examples from the Tarjan-Yannakakis paper
g1 <- graph.formula(A-B:C:I, B-A:C:D, C-A:B:E:H, D-B:E:F,
                    E-C:D:F:H, F-D:E:G, G-F:H, H-C:E:G:I,
                    I-A:H)
maximum.cardinality.search(g1)
is.chordal(g1, fillin=TRUE)

g2 <- graph.formula(A-B:E, B-A:E:F:D, C-E:D:G, D-B:F:E:C:G,
                    E-A:B:C:D:F, F-B:D:E, G-C:D:H:I, H-G:I:J,
                    I-G:H:J, J-H:I)
maximum.cardinality.search(g2)
is.chordal(g2, fillin=TRUE)
}
\keyword{graphs}
