\name{is_separator}
\alias{is.separator}
\alias{is_separator}
\alias{is_min_separator}
\alias{is.minimal.separator}
\concept{Vertex separator}
\title{Vertex separators}
\description{These functions check whether a given set of vertices is a
  vertex separator, or a minimal vertex separator.
}
\usage{
is_separator(graph, candidate)
is_min_separator(graph, candidate)
}
\arguments{
  \item{graph}{The input graph. It may be directed, but edge directions
    are ignored.}
  \item{candidate}{A numeric vector giving the vertex ids of the
    candidate separator.}
}
\details{
  \code{is_separator} decides whether the supplied vertex set is a
  vertex separator. A vertex set is a vertex separator if its removal
  results a disconnected graph.

  \code{is_min_separator} decides whether the supplied vertex set is
  a minimal vertex separator. A minimal vertex separator is a vertex
  separator, such that none of its subsets is a vertex separator.

  In the special case of a fully connected graph with \eqn{n} vertices,
  each set of \eqn{n-1} vertices is considered to be a vertex
  separator.
}
\value{
  A logical scalar, whether the supplied vertex set is a (minimal)
  vertex separator or not.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{
  \code{\link{min_separators}} lists all vertex separator of
  minimum size.
}
\examples{
# The graph from the Moody-White paper
mw <- graph_from_literal(1-2:3:4:5:6, 2-3:4:5:7, 3-4:6:7, 4-5:6:7,
                5-6:7:21, 6-7, 7-8:11:14:19, 8-9:11:14, 9-10,
                10-12:13, 11-12:14, 12-16, 13-16, 14-15, 15-16,
                17-18:19:20, 18-20:21, 19-20:22:23, 20-21,
                21-22:23, 22-23)

# Cohesive subgraphs
mw1 <- induced_subgraph(mw, as.character(c(1:7, 17:23)))
mw2 <- induced_subgraph(mw, as.character(7:16))
mw3 <- induced_subgraph(mw, as.character(17:23))
mw4 <- induced_subgraph(mw, as.character(c(7,8,11,14)))
mw5 <- induced_subgraph(mw, as.character(1:7))

check.sep <- function(G) {
  sep <- min_separators(G)
  sapply(sep, is_min_separator, graph=G)
}

check.sep(mw)
check.sep(mw1)
check.sep(mw2)
check.sep(mw3)
check.sep(mw4)
check.sep(mw5)
}
\keyword{graphs}
