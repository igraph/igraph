
library(igraph)

g1 <- graph.ring(10)
g2 <- graph.star(11, center=11, mode="undirected")
gu <- graph.union(g1, g2)
gu

gdu <- graph.disjoint.union(g1, g2)
gdu

####

graph.difference(gu, g1)

####

graph.intersection(gu, g2)

####

graph.complementer(g2)

####

graph.compose(gu, g1)

