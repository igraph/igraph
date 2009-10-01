\name{leading.eigenvector.community}
\alias{leading.eigenvector.community}
\alias{leading.eigenvector.community.step}
\alias{leading.eigenvector.community.naive}
\alias{community.le.to.membership}
\concept{Community structure}
\title{Community structure detecting based on the leading eigenvector
  of the community matrix}
\description{These functions try to find densely connected subgraphs in
  a graph by calculating the leading non-negative eigenvector of the
  modularity matrix of the graph.}
\usage{
leading.eigenvector.community(graph, steps = -1,
       options = igraph.arpack.default)
leading.eigenvector.community.naive(graph, steps = -1,
       options = igraph.arpack.default)
leading.eigenvector.community.step (graph, fromhere = NULL,
    membership = rep(0, vcount(graph)), community = 0,
    options=igraph.arpack.default)
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
  \item{naive}{Logical constant, it defines how the algorithm tries to
    find more divisions after the first division was made. If
    \code{TRUE} then it simply considers both communities as separate
    graphs and then creates modularity matrices for both communities,
    etc. This method however does not maximize modularity, see the paper
    in the references section below. If it is \code{FALSE} then the
    proper method is used which maximizes modularity.}
  \item{fromhere}{An object returned by a previous call to
    \code{leading.eigenvector.community.step}. This will serve as a
    starting point to take another step. This argument is ignored if it
    is \code{NULL}.}
  \item{membership}{The starting community structure. It is a numeric
    vector defining the membership of every vertex in the graph with a
    number between 0 and the total number of communities at this level
    minus one. By default we start with a single community containing
    all vertices. This argument is ignored if \code{fromhere} is not
    \code{NULL}.

    For \code{community.le.to.memberhip} the starting community
    structure on which \code{steps} merges are performed.
  }
  \item{community}{The id of the community which the algorithm will try
    to split.}
  \item{options}{A named list to override some ARPACK options.}
  \item{merges}{The merge matrix, possible from the result of
    \code{leading.eigenvector.community}.}
}
\details{
  The functions documented in these section implement the 
  \sQuote{leading eigenvector} method developed by Mark Newman and 
  published in MEJ Newman: Finding community structure using the
  eigenvectors of matrices, arXiv:physics/0605087. TODO: proper
  citation.
  
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

  \code{leading.eigenvector.community} is the proper implementation of
  the proposed algorithm. \code{leading.eigenvector.community.naive} is
  considered worse, in this implementation a community found after a
  division is considered as a separate graph for further divisions.
  \code{leading.eigenvector.community.step} is the proper
  implementation, but makes only one step, by trying to split the
  specified community.

  From igraph 0.5 these functions use ARPACK to calculate the
  eigenvectors. See \code{\link{arpack}} for details.

  \code{community.le.to.memberhip} creates a membership vector from the
  result of \code{leading.eigenvector.community} or
  \code{leading.eigenvector.community.naive}. It takes \code{membership}
  and permformes \code{steps} merges, according to the supplied
  \code{merges} matrix.
}
\value{
  \code{leading.eigenvector.community} and
  \code{leading.eigenvector.community.naive} return a named list with
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

  \code{leading.eigenvector.community.step} returns a named list with
  the following members:
  \item{membership}{The new membership vector after the split. If no
    split was done then this is the same as the input membership
    vector.}
  \item{split}{Logical scalar, if \code{TRUE} that means that the
    community was successfully splitted.}
  \item{eigenvector}{The eigenvector of the community matrix, or
    \code{NULL} if the \code{eigenvector} argument was \code{FALSE}.}
  \item{eigenvalue}{The largest positive eigenvalue of the modularity
    matrix.}
  \item{options}{Information about the underlying ARPACK computation,
    see \code{\link{arpack}} for details.
  }

  \code{community.le.to.membership} returns a named list with two
  components:
  \item{membership}{A membership vector, a numerical vector
    indication which vertex belongs to which community. The
    communities are always numbered from zero.}
  \item{csize}{A numeric vector giving the sizes of the communities.}
}
\references{
  MEJ Newman: Finding community structure using the eigenvectors of
  matrices, arXiv:physics/0605087
} 
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{modularity}}, \code{\link{walktrap.community}},
  \code{\link{edge.betweenness.community}},
  \code{\link{fastgreedy.community}},
  \code{as.dendrogram} in package \code{stats}.
}
\examples{
g <- graph.full(5) \%du\% graph.full(5) \%du\% graph.full(5)
g <- add.edges(g, c(0,5, 0,10, 5, 10))
leading.eigenvector.community(g)

lec <- leading.eigenvector.community.step(g)
lec$membership
# Try one more split
leading.eigenvector.community.step(g, fromhere=lec, community=0)
leading.eigenvector.community.step(g, fromhere=lec, community=1)
}
\keyword{graphs}
