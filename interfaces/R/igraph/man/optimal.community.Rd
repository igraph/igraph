\name{optimal.community}
\alias{optimal.community}
\concept{Community structure}
\concept{Modularity}
\title{Optimal community structure}
\description{
  This function calculates the optimal community structure of a graph,
  by maximizing the modularity measure over all possible partitions.}
\usage{
optimal.community(graph)
}
\arguments{
  \item{graph}{The input graph. Edge directions are ignored for directed
    graphs.}
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
  \code{optimal.community} returns a \code{\link{communities}}
  object, please see the \code{\link{communities}} manual page for
  details.
}
\references{
  Ulrik Brandes, Daniel Delling, Marco Gaertler, Robert Gorke, Martin
  Hoefer, Zoran Nikoloski, Dorothea Wagner: On Modularity Clustering,
  \emph{IEEE Transactions on Knowledge and Data Engineering}
  20(2):172-188, 2008.
}
\author{Gabor Csardi \email{Gabor.Csardi@unil.ch}}
\seealso{ \code{\link{communities}} for the documentation of the result,
  \code{\link{modularity}}. See also \code{\link{fastgreedy.community}} 
  for a fast greedy optimizer.
}
\examples{
## Zachary's karate club
g <- graph.famous("Zachary")

## We put everything into a big 'try' block, in case 
## igraph was compiled without GLPK support

try({
  ## The calculation only takes a couple of seconds
  oc <- optimal.community(g)

  ## Double check the result
  print(modularity(oc))
  print(modularity(g, membership(oc)))

  ## Compare to the greedy optimizer
  fc <- fastgreedy.community(g)
  print(modularity(fc))
}, silent=TRUE)
}
\keyword{graphs}
