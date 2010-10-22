
library(igraph)
library(graph, warn.conflicts=FALSE)

g <- erdos.renyi.game(100, 5/100)
N <- igraph.to.graphNEL(g)
g2 <- igraph.from.graphNEL(N)
graph.isomorphic.vf2(g, g2)

## Attributes

V(g)$name <- as.character(vcount(g):1)
E(g)$weight <- sample(1:10, ecount(g), replace=TRUE)
g$name <- "Foobar"

N <- igraph.to.graphNEL(g)
g2 <- igraph.from.graphNEL(N)
graph.isomorphic(g, g2)

all(V(g)$name == V(g2)$name)
A <- get.adjacency(g, attr="weight")
A2 <- get.adjacency(g2, attr="weight")
all(A == A2)
g$name == g2$name

