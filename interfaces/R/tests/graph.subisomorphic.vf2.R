
library(igraph)

set.seed(42)
g1 <- erdos.renyi.game(20,6/20)
g2 <- erdos.renyi.game(20,6/20)
g <- g1 %du% g2

graph.subisomorphic.vf2(g, g1)
graph.subisomorphic.vf2(g, g2)

