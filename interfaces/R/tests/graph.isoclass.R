
library(igraph)

g1 <- graph.isocreate(3, 10)
g2 <- graph.isocreate(3, 11)
graph.isoclass(g1)
graph.isoclass(g2)

g1 <- add.vertices(g1, 3)
graph.isoclass.subgraph(g1, 1:3)
graph.isoclass.subgraph(g1 %du% g2, 1:3)
graph.isoclass.subgraph(g1 %du% g2, 7:9)

