\name{walktrap.community}
\alias{walktrap.community}
\concept{Random walk}
\concept{Community structure}
\title{Community strucure via short random walks}
\description{This function tries to find densely connected subgraphs,
  also called communities in a graph via random walks. The idea is that
  short random walks tend to stay in the same community.
}
\usage{
walktrap.community(graph, weights = E(graph)$weight, steps = 4, merges =
          TRUE, modularity = TRUE, labels = TRUE, membership = TRUE)
}
\arguments{
  \item{graph}{The input graph, edge directions are ignored in directed
    graphs.}
  \item{weights}{The edge weights.}
  \item{steps}{The length of the random walks to perform.}
  \item{merges}{Logical scalar, whether to include the merge matrix in
    the result.}
  \item{modularity}{Logical scalar, whether to include the vector of the
    modularity scores in the result. If the \code{membership} argument
    is true, then it will be always calculated.}
  \item{labels}{Logical scalar, if \code{TRUE} and the graph has a
    vertex attribute called \code{name} then it will be included in the
    result, in the list member called \code{labels}.}
  \item{membership}{Logical scalar, whether to calculate the membership
    vector for the split corresponding to the highest modularity value.}
}
\details{
  This function is the implementation of the Walktrap community
  finding algorithm, see Pascal Pons, Matthieu Latapy: Computing
  communities in large networks using random walks, 
  http://arxiv.org/abs/physics/0512106
}
\value{
  A named list with three members:
  \item{merges}{ The merges performed by the
    algorithm will be stored here. Each merge is a
    row in a two-column matrix and contains the ids of the merged
    communities. Communities are numbered from zero and cluster number 
    smaller than the number of nodes in the network belong to the
    individual vertices as singleton communities. In each step a new
    community is created from two other communities and its id will be 
    one larger than the largest community id so far. This means that 
    before the first merge we have \code{n} communities (the number of
    vertices in the graph) numbered from zero to \code{n-1}. The first
    merge created community \code{n}, the second community \code{n+1}, etc.
  }
  \item{modularity}{ Numeric vector, the modularity score of the
    current community structure after each merge operation. 
  }
  \item{labels}{The labels of the vertices in the graph. The
    \code{name} vertex attribute will be copied here, if exists.
  }
  \item{membership}{ If requested, then the membership vector that
    belongs to the best split, in terms of highest modularity score.
  }
}
\references{Pascal Pons, Matthieu Latapy: Computing
  communities in large networks using random walks, 
  http://arxiv.org/abs/physics/0512106
}
\author{Pascal Pons \email{google@for.it} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the R and igraph interface}
\seealso{\code{\link{modularity}} and 
  \code{\link{fastgreedy.community}},
  \code{\link{spinglass.community}},
  \code{\link{leading.eigenvector.community}},
  \code{\link{edge.betweenness.community}} for other community detection
  methods.
}
\examples{
g <- graph.full(5) \%du\% graph.full(5) \%du\% graph.full(5)
g <- add.edges(g, c(0,5, 0,10, 5, 10))
walktrap.community(g)
}
\keyword{graphs}
