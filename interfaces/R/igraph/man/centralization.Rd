\name{centralization}
\alias{centralization}
\alias{centralize.scores}
\alias{centralization.degree}
\alias{centralization.closeness}
\alias{centralization.betweenness}
\alias{centralization.evcent}
\alias{centralization.degree.tmax}
\alias{centralization.closeness.tmax}
\alias{centralization.betweenness.tmax}
\alias{centralization.evcent.tmax}
\title{Centralization of a graph.tmax}
\description{
  Centralization is a method for creating a graph level centralization
  measure from the centrality scores of the vertices.
}
\usage{
centralize.scores (scores, theoretical.max, normalized = TRUE)

centralization.degree (graph, mode = c("all", "out", "in", "total"),
    loops = TRUE, normalized = TRUE) 
centralization.closeness (graph, mode = c("out", "in", "all", "total"),
    normalized = TRUE)
centralization.betweenness (graph, directed = TRUE, nobigint = TRUE,
    normalized = TRUE)
centralization.evcent (graph, directed = FALSE, scale = TRUE,
    options = igraph.arpack.default, normalized = TRUE)

centralization.degree.tmax (graph = NULL, nodes = 0,
    mode = c("all", "out", "in", "total"), loops = FALSE)
centralization.closeness.tmax (graph = NULL, nodes = 0,
    directed = TRUE)
centralization.betweenness.tmax (graph = NULL, nodes = 0,
    directed = TRUE)
centralization.evcent.tmax (graph = NULL, nodes = 0,
    directed = FALSE, scale = TRUE)
}
\arguments{
  \item{scores}{The vertex level centrality scores.}
  \item{theoretical.max}{Real scalar. The graph level centrality score
    of the most centralized graph with the same number of vertices as
    the graph under study. This is only used if the \code{normalized}
    argument is set to \code{TRUE}.}
  \item{normalized}{Logical scalar. Whether to normalize the graph level
    centrality score by dividing the supplied theoretical maximum.}
  \item{graph}{The input graph. For the \dQuote{tmax} functions it can
    be \code{NULL}, see the details below.}
  \item{mode}{This is the same as the \code{mode} argument of
    \code{degree} and \code{closeness}.}
  \item{loops}{Logical scalar, whether to consider loops edges when
    calculating the degree.}
  \item{directed}{logical scalar, whether to use directed shortest paths
    for calculating betweenness.}
  \item{nobigint}{Logical scalar, whether to use big integers for the
    betweenness calculation. This argument is passed to the
    \code{\link{betweenness}} function.}
  \item{scale}{Whether to rescale the eigenvector centrality scores,
    such that the maximum score is one.}
  \item{nodes}{The number of vertices. This is ignored if the graph is
    given.}
  \item{options}{This is passed to \code{\link{evcent}}, the options for
    the ARPACK eigensolver.}
}
\details{
  Centralization is a general method for calculating a graph-level
  centrality score based on node-level centrality measure. The
  formula for this is

  \deqn{C(G)=\sum_v (\max_w c_w - c_v),}{
    C(G)=sum( max(c(w), w) - c(v),v),}
  
  where \eqn{c_v}{c(v)} is the centrality of vertex \eqn{v}.

  The graph-level centrality score can be normalized by dividing by the
  maximum theoretical score for a graph with the same number of
  vertices, using the same parameters, e.g. directedness, whether we
  consider loop edges, etc.

  For degree, closeness and betweenness the most centralized structure
  is some version of the star graph, in-star, out-star or undirected
  star.

  For eigenvector centrality the most centralized structure is the graph
  with a single edge (and potentially many isolates).

  \code{centralize.scores} using the general centralization formula to
  calculate a graph-level score from vertex-level scores.

  \code{centralization.degree}, \code{centralization.closeness},
  \code{centralization.betweenness} calculate both the vertex-level and
  the graph-level indices.

  \code{centralization.degree.tmax},
  \code{centralization.closeness.tmax},
  \code{centralization.betweenness.tmax} and
  \code{centralization.evcent.tmax} return the theoretical maximum
  scores. They operate in two modes. In the first mode, a graph is given
  and the maximum score is calculated based on that. E.g. the number of
  vertices and directedness is taken from this graph.

  The other way to call these functions is to omit the
  \code{graph} argument, but explicitly specify the rest of the
  arguments.
}
\value{
  For \code{centralize.scores} a real scalar.

  For \code{centralization.degree}, \code{centralization.closeness} and
  \code{centralization.betweenness} a named list with the following
  components:
  \item{res}{The node-level centrality scores.}
  \item{centralization}{The graph level centrality index.}
  \item{theoretical_max}{The maximum theoretical graph level
    centralization score for a graph with the given number of vertices,
    using the same parameters. If the \code{normalized} argument was
    \code{TRUE}, then the result was divided by this number.}

  For \code{centralization.evcent} a named list with the following
  components:
  \item{vector}{The node-level centrality scores.}
  \item{value}{The corresponding eigenvalue.}
  \item{options}{ARPACK options, see the return value of
  \code{\link{evcent}} for details.} 
  \item{centralization}{The graph level centrality index.}
  \item{theoretical_max}{The same as above, the theoretical maximum
    centralization score for a graph with the same number of vertices.}

  For \code{centralization.degree.tmax},
  \code{centralization.closeness.tmax},
  \code{centralization.betweenness.tmax} and
  \code{centralization.evcent.tmax} a real scalar.
}
\references{
  Freeman, L.C.  (1979).  Centrality in Social Networks I:
  Conceptual Clarification. \emph{Social Networks} 1, 215--239.

  Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis:
  Methods and Applications.}  Cambridge University Press.
}
\author{Gabor Csardi \email{Gabor.Csardi@unil.ch}}
%\seealso{}
\examples{
# A BA graph is quite centralized
g <- ba.game(1000, m=4)
centralization.degree(g)$centralization
centralization.closeness(g, mode="all")$centralization
centralization.evcent(g, directed=FALSE)$centralization

# The most centralized graph according to eigenvector centrality
g0 <- graph( c(2,1), n=10, dir=FALSE )
g1 <- graph.star(10, mode="undirected")
centralization.evcent(g0)$centralization
centralization.evcent(g1)$centralization
}
\keyword{graphs}

