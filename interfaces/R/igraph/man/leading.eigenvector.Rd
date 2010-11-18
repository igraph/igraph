\name{leading.eigenvector.community}
\alias{leading.eigenvector.community}
\alias{community.le.to.membership}
\concept{Community structure}
\title{Community structure detecting based on the leading eigenvector
  of the community matrix}
\description{This function tries to find densely connected subgraphs in
  a graph by calculating the leading non-negative eigenvector of the
  modularity matrix of the graph.}
\usage{
leading.eigenvector.community(graph, steps = -1, start = NULL,
       options = igraph.arpack.default)
community.le.to.membership(merges, steps, membership) 
}
\arguments{
  \item{graph}{The input graph. Should be undirected as the method needs
    a symmetric matrix.}
  \item{steps}{The number of steps to take, this is actually the number
    of tries to make a step. It is not a particularly useful parameter.

    For \code{community.le.to.membership} the number of merges to
    produce from the supplied \code{membership} vector.
  }
  \item{start}{\code{NULL}, or a numeric membership vector, giving the
    start configuration of the algorithm.}
  \item{membership}{The starting community
    structure on which \code{steps} merges are performed.
  }
  \item{options}{A named list to override some ARPACK options.}
  \item{merges}{The merge matrix, possible from the result of
    \code{leading.eigenvector.community}.}
}
\details{
  The function documented in these section implements the 
  \sQuote{leading eigenvector} method developed by Mark Newman, see the
  reference below.
  
  The heart of the method is the definition of the modularity matrix,
  \code{B}, which is \code{B=A-P}, \code{A} being the adjacency matrix of
  the (undirected) 
  network, and \code{P} contains the probability that certain edges are
  present according to the \sQuote{configuration model}. In
  other words, a \code{P[i,j]} element of \code{P} is the probability that
  there is an edge between vertices \code{i} and \code{j} in a random
  network in which the degrees of all vertices are the same as in the
  input graph.
  
  The leading eigenvector method works by calculating the eigenvector
  of the modularity matrix for the largest positive eigenvalue and
  then separating vertices into two community based on the sign of
  the corresponding element in the eigenvector. If all elements in
  the eigenvector are of the same sign that means that the network
  has no underlying comuunity structure.
  Check Newman's paper to understand why this is a good method for
  detecting community structure.

  \code{community.le.to.memberhip} creates a membership vector from the
  result of \code{leading.eigenvector.community}. It takes
  \code{membership} and permformes \code{steps} merges, according to the
  supplied \code{merges} matrix.
}
\value{
  \code{leading.eigenvector.community} returns a named list with
  the following members:
  \item{membership}{The membership vector at the end of the algorithm,
    when no more splits are possible.}
  \item{merges}{The merges matrix starting from the state described by
    the \code{membership} member. This is a two-column matrix and each
    line describes a merge of two communities, the first line is the
    first merge and it creates community \sQuote{\code{N}}, \code{N}
    is the number of initial communities in the graph, the second line
    creates community \code{N+1}, etc.
  }
  \item{options}{Information about the underlying ARPACK computation,
    see \code{\link{arpack}} for details.
  }

  \code{community.le.to.membership} returns a named list with two
  components:
  \item{membership}{A membership vector, a numerical vector
    indication which vertex belongs to which community. The
    communities are always numbered from one.}
  \item{csize}{A numeric vector giving the sizes of the communities.}
}
\references{
  MEJ Newman: Finding community structure using the eigenvectors of
  matrices, Physical Review E 74 036104, 2006.
} 
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{modularity}}, \code{\link{walktrap.community}},
  \code{\link{edge.betweenness.community}},
  \code{\link{fastgreedy.community}},
  \code{\link[stats]{as.dendrogram}}
}
\examples{
g <- graph.full(5) \%du\% graph.full(5) \%du\% graph.full(5)
g <- add.edges(g, c(1,6, 1,11, 6, 11))
lec <- leading.eigenvector.community(g)
lec

leading.eigenvector.community(g, start=membership(lec))
}
\keyword{graphs}
