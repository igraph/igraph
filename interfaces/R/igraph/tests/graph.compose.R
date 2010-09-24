
library(igraph)

g1 <- erdos.renyi.game(50, 3/50, directed=TRUE)
gi <- graph( rep(1:vcount(g1), each=2), dir=TRUE )

g2 <- graph.compose(g1, gi)
g3 <- graph.compose(gi, g1)

graph.isomorphic(g1, g2)
graph.isomorphic(g1, g3)

############

