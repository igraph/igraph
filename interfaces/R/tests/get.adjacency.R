
library(igraph)

g <- erdos.renyi.game(50, 1/50)
A <- get.adjacency(g, sparse=FALSE)
g2 <- graph.adjacency(A, mode="undirected")
graph.isomorphic(g, g2)

###

A <- get.adjacency(g, sparse=TRUE)
g2 <- graph.adjacency(A, mode="undirected")
graph.isomorphic(g, g2)

###

g <- erdos.renyi.game(50, 2/50, directed=TRUE)
A <- get.adjacency(g, sparse=FALSE)
g2 <- graph.adjacency(A)
graph.isomorphic(g, g2)

###

A <- get.adjacency(g, sparse=TRUE)
g2 <- graph.adjacency(A)
graph.isomorphic(g, g2)

