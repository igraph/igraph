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
  \item{graph}{The input graph.}
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
  \code{walktrap.community} returns a \code{\link{communities}}
  object, please see the \code{\link{communities}} manual page for
  details.
}
\references{Pascal Pons, Matthieu Latapy: Computing
  communities in large networks using random walks, 
  http://arxiv.org/abs/physics/0512106
}
\author{Pascal Pons \email{google@for.it} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the R and igraph interface}
\seealso{
  See \code{\link{communities}} on getting the actual membership vector,
  merge matrix, modularity score, etc.
  
  \code{\link{modularity}} and 
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
