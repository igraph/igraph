\name{graph.maxflow}
\alias{graph.maxflow}
\title{Maximum flow in a network}
\description{In a graph where each edge has a given flow capacity the
  maximal flow between two vertices is calculated.}
\usage{
graph.maxflow(graph, source, target, capacity=NULL)
}
\arguments{
  \item{graph}{The input graph.}
  \item{source}{The id of the source vertex.}
  \item{target}{The id of the target vertex (sometimes also called sink).}
  \item{capacity}{Vector giving the capacity of the edges. If this is
    \code{NULL} (the default) then the \code{capacity} edge attribute is
    used.}
}
% \details{}
\value{A numeric constant, the maximum flow.}
\references{
  A. V. Goldberg and R. E. Tarjan: ``A New Approach to the Maximum Flow
  Problem'' Journal of the ACM 35:921-940, 1988.  
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{shortest.paths}}}
\examples{
}
\keyword{graphs}
