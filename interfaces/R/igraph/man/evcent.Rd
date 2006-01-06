\name{evcent}
\alias{evcent}
\title{Find Eigenvector Centrality Scores of Network Positions}
\description{
  \code{evcent} takes a graph (\code{graph}) and returns the
  eigenvector centralities of positions \code{v} within it
}
\usage{
evcent(graph, v=c(1:vcount(graph))
}
\arguments{
  \item{graph}{Graph to be analyzed.}
  \item{v}{Numeric vector indicating which vertices are to be included
           in the calculation.  By default, all vertices are included. }
}
\details{
Eigenvector centrality scores correspond to the values of the first
eigenvector of the graph adjacency matrix; these scores may, in turn, be
interpreted as arising from a reciprocal process in which the centrality
of each actor is proportional to the sum of the centralities of those
actors to whom he or she is connected.  In general, vertices with high
eigenvector centralities are those which are connected to many other
vertices which are, in turn, connected to many others (and so on).  (The
perceptive may realize that this implies that the largest values will be
obtained by individuals in large cliques (or high-density
substructures).  This is also intelligible from an algebraic point of
view, with the first eigenvector being closely related to the best
rank-1 approximation of the adjacency matrix (a relationship which is
easy to see in the special case of a diagonalizable symmetric real
matrix via the \eqn{SLS^-1}{$S \Lambda S^{-1}$} decomposition).)  
}
\value{
  A vector containing the centrality scores.
}
\references{ Bonacich, P.  (1987).  ``Power and Centrality: A Family of
Measures.'' \emph{American Journal of Sociology}, 92, 1170-1182. 

Katz, L.  (1953).  ``A New Status Index Derived from Sociometric
Analysis.''  \emph{Psychometrika}, 18, 39-43. 
 }
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu} and Carter T. Butts
\email{buttsc@uci.edu} for the manual page} 
\section{WARNING }{\code{evcent} will not symmetrize your data before
extracting eigenvectors; don't send this routine asymmetric matrices
unless you really mean to do so.}
% \seealso{}
\examples{
#Generate some test data
g <- graph.ring(10, directed=FALSE)
#Compute eigenvector centrality scores
evcent(g)
}
\keyword{graphs}
