\name{decompose.graph}
\alias{decompose.graph}
\concept{Graph decomposition}
\concept{Graph component}
\title{Decompose a graph into components}
\description{Creates a separate graph for each component of a graph.}
\usage{
decompose.graph(graph, mode = c("weak", "strong"),
       max.comps = NA, min.vertices = 0) 
}
\arguments{
  \item{graph}{The original graph.}
  \item{mode}{Character constant giving the type of the components,
    wither \code{weak} for weakly connected components or \code{strong}
    for strongly connected components.}
  \item{max.comps}{The maximum number of components to return. The first
    \code{max.comps} components will be returned (which hold at least
    \code{min.vertices} vertices, see the next parameter), the others
    will be ignored. Supply \code{NA} here if you don't want to limit
    the number of components.}
  \item{min.vertices}{The minimum number of vertices a component should
    contain in order to place it in the result list. Eg. supply 2
    here to ignore isolate vertices.}
}
% \details{}
\value{
  A list of graph objects.
}
%\references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{is.connected}} to decide whether a graph is
  connected, \code{\link{clusters}} to calculate the connected
  components of a graph.}
\examples{
# the diameter of each component in a random graph
g <- erdos.renyi.game(1000, 1/1000)
comps <- decompose.graph(g, min.vertices=2)
sapply(comps, diameter)
}
\keyword{graphs}
