
library(igraph)

g <- erdos.renyi.game(100, 3/100)
edges <- unlist(lapply(seq_len(ecount(g)), get.edge, graph=g))
g2 <- graph(edges, dir=FALSE)
graph.isomorphic(g, g2)

