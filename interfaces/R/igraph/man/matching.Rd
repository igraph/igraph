\name{graph.matching}
\alias{is.matching}
\alias{is.maximal.matching}
\alias{maximum.bipartite.matching}
\concept{Matching}
\concept{Maximal matching}
\concept{Maximum bipartite matching}
\title{Graph matching}
\description{A matching in a graph means the selection of a set of edges
  that are pairwise non-adjacenct, i.e. they have no common incident
  vertices. A matching is maximal if it is not a proper subset of any
  other matching.}
\usage{
is.matching(graph, matching, types = NULL)
is.maximal.matching(graph, matching, types = NULL)
maximum.bipartite.matching(graph, types = NULL, weights = NULL,
       eps = .Machine$double.eps)
}
\arguments{
  \item{graph}{The input graph. It might be directed, but edge directions
    will be ignored.}
  \item{types}{Vertex types, if the graph is bipartite. By default they
    are taken from the \sQuote{\code{type}} vertex attribute, if
    present.}
  \item{matching}{A potential matching. An integer vector that gives the
    pair in the matching for each vertex. For vertices without a pair,
    supply \code{NA} here.}
  \item{weights}{Potential edge weights. If the graph has an edge
    attribute called \sQuote{\code{weight}}, and this argument is
    \code{NULL}, then the edge attribute is used automatically.}
  \item{eps}{A small real number used in equality tests in the weighted
    bipartite matching algorithm. Two real numbers are considered
    equal in the algorithm if their difference is smaller than
    \code{eps}. This is required to avoid the accumulation of numerical
    errors. By default it is set to the smallest \eqn{x}, such that
    \eqn{1+x \ne 1}{1+x != 1} holds. If you are running the algorithm
    with no  weights, this argument is ignored.}
}
\details{
  \code{is.matching} checks a matching vector and verifies whether its
  length matches the number of vertices in the given graph, its values
  are between zero (inclusive) and the number of vertices (inclusive),
  and whether there exists a corresponding edge in the graph for every
  matched vertex pair. For bipartite graphs, it also verifies whether
  the matched vertices are in different parts of the graph.

  \code{is.maximal.matching} checks whether a matching is maximal.
  A matching is maximal if and only if there exists no unmatched vertex
  in a graph such that one of its neighbors is also unmatched.

  \code{maximum.bipartite.matching} calculates a maximum matching in a
  bipartite graph. A matching in a bipartite graph is a partial
  assignment of vertices of the first kind to vertices of the second
  kind such that each vertex of the first kind is matched to at most one
  vertex of the second kind and vice versa, and matched vertices must be
  connected by an edge in the graph. The size (or cardinality) of a
  matching is the number of edges. A matching is a maximum matching if
  there exists no other matching with larger cardinality. For weighted
  graphs, a maximum matching is a matching whose edges have the largest
  possible total weight among all possible matchings.
 
  Maximum matchings in bipartite graphs are found by the push-relabel
  algorithm with greedy initialization and a global relabeling after
  every \eqn{n/2} steps where \eqn{n} is the number of vertices in the
  graph.
}
\value{
  \code{is.matching} and \code{is.maximal.matching} return a logical
  scalar.

  \code{maximum.bipartite.matching} returns a list with components:
  \item{matching_size}{The size of the matching, i.e. the number of
    edges connecting the matched vertices.}
  \item{matching_weight}{The weights of the matching, if the graph was
    weighted. For unweighted graphs this is the same as the size of the
    matching.}
  \item{matching}{The matching itself. Numeric vertex id, or vertex
    names if the graph was named. Non-matched vertices are denoted by
    \code{NA}.}
}
% \references{}
\author{Tamas Nepusz \email{ntamas@gmail.com}}
\examples{
g <- graph.formula( a-b-c-d-e-f )
m1 <- c("b", "a", "d", "c", "f", "e")   # maximal matching
m2 <- c("b", "a", "d", "c", NA, NA)     # non-maximal matching
m3 <- c("b", "c", "d", "c", NA, NA)     # not a matching
is.matching(g, m1)
is.matching(g, m2)
is.matching(g, m3)
is.maximal.matching(g, m1)
is.maximal.matching(g, m2)
is.maximal.matching(g, m3)

V(g)$type <- c(FALSE,TRUE)
str(g, v=TRUE)
maximum.bipartite.matching(g)

g2 <- graph.formula( a-b-c-d-e-f-g )
V(g2)$type <- rep(c(FALSE,TRUE), length=vcount(g2))
str(g2, v=TRUE)
maximum.bipartite.matching(g2)
}
\keyword{graphs}
