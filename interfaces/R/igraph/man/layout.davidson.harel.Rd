\name{layout.davidson.harel}
\alias{layout.davidson.harel}
\concept{Graph layout}
\title{The Davidson-Harel layout algorithm}
\description{Place vertices of a graph on the plane, according to the
  simulated annealing algorithm by Davidson and Harel.}
\usage{
layout.davidson.harel (graph, coords=NULL, maxiter=10,
           fineiter=max(10, log2(vcount(graph))), cool.fact=0.75,
           weight.node.dist=1.0, weight.border=0.0,
           weight.edge.lengths=graph.density(graph) / 10,
           weight.edge.crossings=1.0 - sqrt(graph.density(graph)),
           weight.node.edge.dist=0.2 * (1-graph.density(graph)))
}
\arguments{
  \item{graph}{The graph to lay out. Edge directions are ignored.}
  \item{coords}{Optional starting positions for the vertices. If this
    argument is not \code{NULL} then it should be an appropriate
    matrix of starting coordinates.}
  \item{maxiter}{Number of iterations to perform in the first phase.}
  \item{fineiter}{Number of iterations in the fine tuning phase.}
  \item{cool.fact}{Cooling factor.}
  \item{weight.node.dist}{Weight for the node-node distances
    component of the energy function. }
  \item{weight.border}{Weight for the distance from the border
    component of the energy function. It can be set to zero, if
    vertices are allowed to sit on the border.}
  \item{weight.edge.lengths}{Weight for the edge length component
    of the energy function.}
  \item{weight.edge.crossings}{Weight for the edge crossing component
    of the energy function.}
  \item{weight.node.edge.dist}{Weight for the node-edge distance
    component of the energy function.}
}
\details{
  This function implements the algorithm by Davidson and Harel,
  see Ron Davidson, David Harel: Drawing Graphs Nicely Using 
  Simulated Annealing. ACM Transactions on Graphics 15(4),
  pp. 301-331, 1996.
  
  The algorithm uses simulated annealing and a sophisticated 
  energy function, which is unfortunately hard to parameterize
  for different graphs. The original publication did not disclose any
  parameter values, and the ones below were determined by 
  experimentation.
  
  The algorithm consists of two phases, an annealing phase, and a 
  fine-tuning phase. There is no simulated annealing in the second
  phase.
  
  Our implementation tries to follow the original publication, as
  much as possible. The only major difference is that coordinates are
  explicitly kept within the bounds of the rectangle of the layout.
}
\value{
  A two- or three-column matrix, each row giving the coordinates of a
  vertex, according to the ids of the vertex ids.
}
\references{
  Ron Davidson, David Harel: Drawing Graphs Nicely Using 
  Simulated Annealing. \emph{ACM Transactions on Graphics} 15(4),
  pp. 301-331, 1996.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout.fruchterman.reingold}},
  \code{\link{layout.kamada.kawai}} for other layout algorithms.}
\examples{
set.seed(42)
L <- layout.davidson.harel
## Figures from the paper
g_1b <- graph.star(19, mode="undirected") + path(c(2:19, 2)) +
  path(c(seq(2, 18, by=2), 2))
plot(g_1b, layout=L)

g_2 <- graph.lattice(c(8, 3)) + edges(1,8, 9,16, 17,24)
plot(g_2, layout=L)

g_3 <- graph.empty(n=70)
plot(g_3, layout=L)

g_4 <- graph.empty(n=70, directed=FALSE) + edges(1:70)
plot(g_4, layout=L, vertex.size=5, vertex.label=NA)

g_5a <- graph.ring(24)
plot(g_5a, layout=L, vertex.size=5, vertex.label=NA)

g_5b <- graph.ring(40)
plot(g_5b, layout=L, vertex.size=5, vertex.label=NA)

g_6 <- graph.lattice(c(2,2,2))
plot(g_6, layout=L)

g_7 <- graph.formula(1:3:5 -- 2:4:6)
plot(g_7, layout=L, vertex.label=V(g_7)$name)

g_8 <- graph.ring(5) + graph.ring(10) + graph.ring(5) +
  edges(1,6, 2,8, 3, 10, 4,12, 5,14,
        7,16, 9,17, 11,18, 13,19, 15,20)
plot(g_8, layout=L, vertex.size=5, vertex.label=NA)

g_9 <- graph.lattice(c(3,2,2))
plot(g_9, layout=L, vertex.size=5, vertex.label=NA)

g_10 <- graph.lattice(c(6,6))
plot(g_10, layout=L, vertex.size=5, vertex.label=NA)

g_11a <- graph.tree(31, 2, mode="undirected")
plot(g_11a, layout=L, vertex.size=5, vertex.label=NA)

g_11b <- graph.tree(21, 4, mode="undirected")
plot(g_11b, layout=L, vertex.size=5, vertex.label=NA)

g_12 <- graph.empty(n=37, directed=FALSE) +
  path(1:5,10,22,31,37:33,27,16,6,1) + path(6,7,11,9,10) + path(16:22) +
  path(27:31) + path(2,7,18,28,34) + path(3,8,11,19,29,32,35) +
  path(4,9,20,30,36) + path(1,7,12,14,19,24,26,30,37) +
  path(5,9,13,15,19,23,25,28,33) + path(3,12,16,25,35,26,22,13,3)
plot(g_12,  layout=L, vertex.size=5, vertex.label=NA)
}
