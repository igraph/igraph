\name{canonical.permutation}
\alias{canonical.permutation}
\concept{Canonical permutation}
\concept{BLISS}
\title{Canonical permutation of a graph}
\description{The canonical permutation brings every isomorphic graphs
  into the same (labeled) graphs.}
\usage{
canonical.permutation(graph, sh="fm")
}
\arguments{
  \item{graph}{The input graph, treated as undirected.}
  \item{sh}{Type of the heuristics to use for the BLISS
    algorithm. See details for possible values.}
}
\details{
  \code{canonical.permutation} computes a permutation which brings the
  graph into canonical form, as defined by the BLISS algorithm.
  All isomorphic graphs have the same canonical form.

  See the paper below for the details about BLISS. This and more
  information is available at
  \url{http://www.tcs.hut.fi/Software/bliss/index.html}.

  The possible values for the \code{sh} argument are:
  \itemize{
    \item{f}{First non-singleton cell.}
    \item{fl}{First largest non-singleton cell.}
    \item{fs}{First smallest non-singleton cell.}
    \item{fm}{First maximally non-trivially connectec non-singleton
      cell.}
    \item{flm}{Largest maximally non-trivially connected non-singleton
      cell.}
    \item{fsm}{Smallest maximally non-trivially connected non-singleton
      cell.}
  }
  See the paper in references for details about these.
}
\value{
  A list with the following members:
  \item{labeling}{The canonical parmutation which takes the input
    graph into canonical form. A numeric vector, the first element is
    the new label of vertex 0, the second element for vertex 1, etc. }
  \item{info}{Some information about the BLISS computation. A named
    list with the following members:
    \itemize{
      \item{nof\_nodes}{The number of nodes in the search tree.}
      \item{nof\_leaf\_nodes}{The number of leaf nodes in the search
	tree.}
      \item{nof\_bad\_nodes}{Number of bad nodes.}
      \item{nof\_canupdates}{Number of canrep updates.}
      \item{max\_level}{Maximum level.}
      \item{group\_size}{The size of the automorphism group of the
	input graph, as a string. This number is exact if igraph was
	compiled with the GMP library, and approximate otherwise.}
    }
  }
}
\references{
  Tommi Junttila and Petteri Kaski: Engineering an Efficient Canonical
  Labeling Tool for Large and Sparse Graphs, \emph{Proceedings of the
    Ninth Workshop on Algorithm Engineering and Experiments and the
    Fourth Workshop on Analytic Algorithms and Combinatorics.} 2007.
}
\author{Tommi Junttila for BLISS, Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the igraph and R interfaces.}
\seealso{\code{\link{permute.vertices}} to apply a permutation to a
  graph, \code{\link{graph.isomorphic}} for deciding graph isomorphism,
  possibly based on canonical labels.}
\examples{
## Calculate the canonical form of a random graph
g1 <- erdos.renyi.game(10, 20, type="gnm")
cp1 <- canonical.permutation(g1)
cf1 <- permute.vertices(g1, cp1$labeling)

## Do the same with a random permutation of it
g2 <- permute.vertices(g1, sample(vcount(g1))-1)
cp2 <- canonical.permutation(g2)
cf2 <- permute.vertices(g2, cp2$labeling)

## Check that they are the same
el1 <- get.edgelist(cf1)
el2 <- get.edgelist(cf2)
el1 <- el1[ order(el1[,1], el1[,2]), ]
el2 <- el2[ order(el2[,1], el2[,2]), ]
all(el1 == el2)
}
\keyword{graphs}
