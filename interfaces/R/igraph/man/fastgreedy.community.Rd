\name{fastgreedy.community}
\alias{fastgreedy.community}
\concept{Community structure}
\concept{Fast greedy community detection}
\title{Community structure via greedy optimization of modularity}
\description{
  This function tries to find dense subgraph, also called communities in
  graphs via directly optimizing a modularity score.
}
\usage{
fastgreedy.community(graph, merges=TRUE, modularity=TRUE,
      weights=E(graph)$weight)
}
\arguments{
  \item{graph}{The input graph}
  \item{merges}{Logical scalar, whether to return the merge matrix.}
  \item{modularity}{Logcal scalar, whether to return a vector containing
    the modularity after each merge.}
  \item{weights}{If not \code{NULL}, then a numeric vector of edge
    weights. The length must match the number of edges in the graph.
    By default the \sQuote{\code{weight}} edge attribute is used as
    weights. If it is not present, then all edges are considered to have
    the same weight.
  }
}
\details{
 This function implements the fast greedy modularity optimization
 algorithm for finding community structure, see 
 A Clauset, MEJ Newman, C Moore: Finding community structure in very
 large networks, http://www.arxiv.org/abs/cond-mat/0408187 for the
 details.
}
\value{
  A named list with the following members:
  \item{merges}{A matrix with two column, this represents a dendogram
    and contains all the merges the algorithm performed. Each line is
    one merge and it is given by the ids of the two communities
    merged. The community ids are integer numbers starting from zero
    and the communities between zero and the number of vertices
    (\code{N}) minus one belong to individual vertices. The first line
    of the matrix gives the first merge, this merge creates community
    \code{N}, the number of vertices, the second merge creates
    community \code{N+1}, etc.
  }
  \item{modularity}{A numeric vector containing the modularity value 
    of the community structure after performing every merge.
  }
}
\references{
  A Clauset, MEJ Newman, C Moore: Finding community structure in very
  large networks, http://www.arxiv.org/abs/cond-mat/0408187
}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the R interface.
}
\seealso{ \code{\link{walktrap.community}},
  \code{\link{spinglass.community}},
  \code{\link{leading.eigenvector.community}},
  \code{\link{edge.betweenness.community}}
}
\examples{
g <- graph.full(5) \%du\% graph.full(5) \%du\% graph.full(5)
g <- add.edges(g, c(0,5, 0,10, 5, 10))
fastgreedy.community(g)
# The highest value of modularity is before performing the last two
# merges. So this network naturally has three communities.
}
\keyword{graphs}
