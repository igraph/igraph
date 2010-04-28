\name{optimal.community}
\alias{optimal.community}
\concept{Community structure}
\concept{Modularity}
\title{Optimal community structure}
\description{
  This function calculates the optimal community structure of a graph,
  by maximizing the modularity measure over all possible partitions.}
\usage{
optimal.community(graph, verbose = igraph.par("verbose"))
}
\arguments{
  \item{graph}{The input graph. Edge directions are ignored for directed
    graphs.}
  \item{verbose}{Logical scalar, whether to print status messages about
    the progress of the algorithm to the screen. Currently not
    implemented.}
}
\details{
  This function calculates the optimal community structure for a
  graph, in terms of maximal modularity score. 
  
  The calculation is done by transforming the modularity maximization
  into an integer programming problem, and then calling the GLPK
  library to solve that. Please the reference below for details. 

  Note that modularity optimization is an NP-complete problem, and
  all known algorithms for it have exponential time complexity. This
  means that you probably don't want to run this function on larger
  graphs. Graphs with up to fifty vertices should be fine, graphs
  with a couple of hundred vertices might be possible.
}
\value{
  A named list with entries:
  \item{modularity}{A numeric scalar, the maximal modularity score for
    the graph.
  }
  \item{membership}{A numeric vector containing integers, one for each
    vertex in the input graph. Vertices labeled with the same integer
    belong to the same community. The community labels start from zero.
  }
}
\references{
  Ulrik Brandes, Daniel Delling, Marco Gaertler, Robert Gorke, Martin
  Hoefer, Zoran Nikoloski, Dorothea Wagner: On Modularity Clustering,
  \emph{IEEE Transactions on Knowledge and Data Engineering}
  20(2):172-188, 2008.
}
\author{Gabor Csardi \email{Gabor.Csardi@unil.ch}}
\seealso{ \code{\link{modularity}}, \code{\link{fastgreedy.community}}
  for a fast greedy optimizer.
}
\examples{
## Zachary's karate club
e <- c( 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8,
       0,10, 0,11, 0,12, 0,13, 0,17, 0,19, 0,21, 0,31,
       1, 2, 1, 3, 1, 7, 1,13, 1,17, 1,19, 1,21, 1,30,
       2, 3, 2, 7, 2,27, 2,28, 2,32, 2, 9, 2, 8, 2,13,
       3, 7, 3,12, 3,13, 4, 6, 4,10, 5, 6, 5,10, 5,16,
       6,16, 8,30, 8,32, 8,33, 9,33,13,33,14,32,14,33,
       15,32,15,33,18,32,18,33,19,33,20,32,20,33,
       22,32,22,33,23,25,23,27,23,32,23,33,23,29,
       24,25,24,27,24,31,25,31,26,29,26,33,27,33,
       28,31,28,33,29,32,29,33,30,32,30,33,31,32,31,33,
       32,33)
g <- graph(e, dir=FALSE)

## We put everything into a big 'try' block, in case 
## igraph was compiled without GLPK support

try({
  ## The calculation only takes a couple of seconds
  oc <- optimal.community(g)

  ## Double check the result
  oc$modularity
  modularity(g, oc$membership)

  ## Compare to the greedy optimizer
  fc <- fastgreedy.community(g)
  max(fc$modularity)
}, silent=TRUE)
}
\keyword{graphs}
