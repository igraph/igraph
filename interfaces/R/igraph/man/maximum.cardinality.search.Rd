\name{maximum.cardinality.search}
\alias{maximum.cardinality.search}
\concept{maximum cardinality search}
\concept{graph decomposition}
\concept{chordal graph}
\title{Maximum cardinality search}
\description{Maximum cardinality search is a simple ordering a vertices
  that is useful in determining the chordality of a graph.}
\usage{
maximum.cardinality.search(graph)
}
\arguments{
  \item{graph}{The input graph. It may be directed, but edge directions
    are ignored, as the algorithm is defined for undirected graphs.}
}
\details{
  Maximum cardinality search visits the vertices in such an order that
  every time the vertex with the most already visited neighbors is
  visited. Ties are broken randomly.

  The algorithm provides a simple basis for deciding whether a graph is
  chordal, see References below, and also \code{\link{is.chordal}}.
}
\value{A list with two components:
  \item{alpha}{Numeric vector. The vertices ordered according to the maximum cardinality
    search.}
  \item{alpham1}{Numeric vector. The inverse of \code{alpha}.}
}
\references{ Robert E Tarjan and Mihalis Yannakakis. (1984). Simple
  linear-time algorithms to test chordality of graphs, test acyclicity
  of hypergraphs, and selectively reduce acyclic hypergraphs.
  \emph{SIAM Journal of Computation} 13, 566--579.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{is.chordal}} }
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
