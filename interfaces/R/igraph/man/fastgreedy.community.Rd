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
      membership=TRUE, weights=E(graph)$weight)
}
\arguments{
  \item{graph}{The input graph}
  \item{merges}{Logical scalar, whether to return the merge matrix.}
  \item{modularity}{Logical scalar, whether to return a vector containing
    the modularity after each merge.}
  \item{membership}{Logical scalar, whether to calculate the membership
    vector corresponding to the maximum modularity score, considering
    all possible community structures along the merges.}
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
  \code{fastgreedy.community} returns a \code{\link{communities}}
  object, please see the \code{\link{communities}} manual page for
  details.
}
\references{
  A Clauset, MEJ Newman, C Moore: Finding community structure in very
  large networks, http://www.arxiv.org/abs/cond-mat/0408187
}
\author{Tamas Nepusz \email{ntamas@gmail.com} and Gabor Csardi
  \email{csardi.gabor@gmail.com} for the R interface.
}
\seealso{
  \code{\link{communities}} for extracting the results.
  
  See also \code{\link{walktrap.community}},
  \code{\link{spinglass.community}},
  \code{\link{leading.eigenvector.community}} and 
  \code{\link{edge.betweenness.community}} for other methods.
}
\examples{
g <- graph.full(5) \%du\% graph.full(5) \%du\% graph.full(5)
g <- add.edges(g, c(1,6, 1,11, 6, 11))
fc <- fastgreedy.community(g)
membership(fc)
sizes(fc)
}
\keyword{graphs}
