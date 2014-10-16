\name{centralization}
\alias{centralization}
\alias{centralize.scores}
\alias{centralize}
\alias{centralization.degree}
\alias{centr_degree}
\alias{centralization.closeness}
\alias{centr_clo}
\alias{centralization.betweenness}
\alias{centr_betw}
\alias{centralization.evcent}
\alias{centr_eigen}
\alias{centralization.degree.tmax}
\alias{centr_degree_tmax}
\alias{centralization.closeness.tmax}
\alias{centr_clo_tmax}
\alias{centralization.betweenness.tmax}
\alias{centr_betw_tmax}
\alias{centralization.evcent.tmax}
\alias{centr_eigen_tmax}
\title{Centralization of a graph.tmax}
\description{
  Centralization is a method for creating a graph level centralization
  measure from the centrality scores of the vertices.
}
\usage{
centralize (scores, theoretical.max, normalized = TRUE)

centr_degree (graph, mode = c("all", "out", "in", "total"),
    loops = TRUE, normalized = TRUE) 
centr_clo (graph, mode = c("out", "in", "all", "total"),
    normalized = TRUE)
centr_betw (graph, directed = TRUE, nobigint = TRUE,
    normalized = TRUE)
centr_eigen (graph, directed = FALSE, scale = TRUE,
    options = arpack_defaults, normalized = TRUE)

centr_degree_tmax (graph = NULL, nodes = 0,
    mode = c("all", "out", "in", "total"), loops = FALSE)
centr_clo_tmax (graph = NULL, nodes = 0,
    mode = c("out", "in", "all", "total"))
centr_betw_tmax (graph = NULL, nodes = 0,
    directed = TRUE)
centr_eigen_tmax (graph = NULL, nodes = 0,
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
  \item{options}{This is passed to \code{\link{eigen_centrality}}, the options for
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

  \code{centralize} using the general centralization formula to
  calculate a graph-level score from vertex-level scores.

  \code{centr_degree}, \code{centr_clo},
  \code{centr_betw} calculate both the vertex-level and
  the graph-level indices.

  \code{centr_degree_tmax},
  \code{centr_clo_tmax},
  \code{centr_betw_tmax} and
  \code{centr_eigen_tmax} return the theoretical maximum
  scores. They operate in two modes. In the first mode, a graph is given
  and the maximum score is calculated based on that. E.g. the number of
  vertices and directedness is taken from this graph.

  The other way to call these functions is to omit the
  \code{graph} argument, but explicitly specify the rest of the
  arguments.
}
\value{
  For \code{centralize} a real scalar.

  For \code{centr_degree}, \code{centr_clo} and
  \code{centr_betw} a named list with the following
  components:
  \item{res}{The node-level centrality scores.}
  \item{centralization}{The graph level centrality index.}
  \item{theoretical_max}{The maximum theoretical graph level
    centralization score for a graph with the given number of vertices,
    using the same parameters. If the \code{normalized} argument was
    \code{TRUE}, then the result was divided by this number.}

  For \code{centr_eigen} a named list with the following
  components:
  \item{vector}{The node-level centrality scores.}
  \item{value}{The corresponding eigenvalue.}
  \item{options}{ARPACK options, see the return value of
  \code{\link{eigen_centrality}} for details.} 
  \item{centralization}{The graph level centrality index.}
  \item{theoretical_max}{The same as above, the theoretical maximum
    centralization score for a graph with the same number of vertices.}

  For \code{centr_degree_tmax},
  \code{centr_clo_tmax},
  \code{centr_betw_tmax} and
  \code{centr_eigen_tmax} a real scalar.
}
\references{
  Freeman, L.C.  (1979).  Centrality in Social Networks I:
  Conceptual Clarification. \emph{Social Networks} 1, 215--239.

  Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis:
  Methods and Applications.}  Cambridge University Press.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
%\seealso{}
\examples{
# A BA graph is quite centralized
g <- sample_pa(1000, m=4)
centr_degree(g)$centralization
centr_clo(g, mode="all")$centralization
centr_eigen(g, directed=FALSE)$centralization

# The most centralized graph according to eigenvector centrality
g0 <- graph( c(2,1), n=10, dir=FALSE )
g1 <- make_star(10, mode="undirected")
centr_eigen(g0)$centralization
centr_eigen(g1)$centralization
}
\keyword{graphs}

