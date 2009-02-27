\name{alpha.centrality}
\alias{alpha.centrality}
\concept{Alpha centrality}
\title{Find Bonacich alpha centrality scores of network positions}
\description{\code{alpha.centrality} calculates the alpha centrality of
  some (or all) vertices in a graph.
}
\usage{
alpha.centrality(graph, nodes=V(graph), alpha=1, loops=FALSE,
                 exo=1, weights=NULL, tol=1e-7)
}
\arguments{
  \item{graph}{The input graph, can be directed or undirected}
  \item{nodes}{Vertex sequence, the vertices for which the alpha
    centrality values are returned. (For technical reasons they
    will be calculated for all vertices first, anyway.)}
  \item{alpha}{Parameter specifying the relative importance of
    endogenous versus exogenous factors in the determination of
    centrality. See details below.}
  \item{loops}{Whether to eliminate loop edges from the graph before the
    calculation.}
  \item{exo}{The exogenous factors, in most cases this is either a
    constant -- the same factor for every node, or a vector giving the
    factor for every vertex. Note that long vectors will be truncated and
    short vectors will be replicated.}
  \item{weights}{Optional positive weight vector for calculating
    weighted closeness. If the graph has a \code{weight} edge
    attribute, then this is used by default.}
  \item{tol}{Tolerance for near-singularities during matrix inversion,
    see \code{\link{solve}}.}
}
\details{
  The alpha centrality measure can be considered as a generalization of
  eigenvector centerality to directed graphs. It was proposed by
  Bonacich in 2001 (see reference below).

  The alpha centrality of the vertices in a graph is defined as the
  solution of the following matrix equation:
  \deqn{x=\alpha A^T x+e,}{x=alpha t(A)x+e,}
  where \eqn{A}{A} is the (not neccessarily symmetric) adjacency matrix of
  the graph, \eqn{e}{e} is the vector of exogenous sources of status of the
  vertices and \eqn{\alpha}{alpha} is the relative importance of the endogenous
  versus exogenous factors.
}
\value{
  A numeric vector contaning the centrality scores for the selected
  vertices.
}
\references{
  Bonacich, P. and Paulette, L. (2001). ``Eigenvector-like measures of
  centrality for asymmetric relations'' \emph{Social Networks}, 23,
  191-201.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}.}
\section{Warning}{Singular adjacency matrices cause problems for this
  algorithm, the routine may fail is certain cases.}
\seealso{\code{\link{evcent}} and \code{\link{bonpow}}}
\examples{
# The examples from Bonacich's paper
g.1 <- graph( c(1,3,2,3,3,4,4,5)-1 )
g.2 <- graph( c(2,1,3,1,4,1,5,1)-1 )
g.3 <- graph( c(1,2,2,3,3,4,4,1,5,1)-1 )
alpha.centrality(g.1)
alpha.centrality(g.2)
alpha.centrality(g.3,alpha=0.5)
}
\keyword{graphs}
