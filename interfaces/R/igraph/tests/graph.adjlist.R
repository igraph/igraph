
library(igraph)

g <- erdos.renyi.game(100, 3/100)
al <- get.adjlist(g)
g2 <- graph.adjlist(al, mode="all")
graph.isomorphic(g, g2)

##

g <- erdos.renyi.game(100, 3/100, dir=TRUE)
al <- get.adjlist(g, mode="out")
g2 <- graph.adjlist(al, mode="out")
graph.isomorphic(g, g2)

