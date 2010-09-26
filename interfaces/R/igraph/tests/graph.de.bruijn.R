
library(igraph)

g <- graph.de.bruijn(2,1)
g2 <- graph.de.bruijn(2,2)
g3 <- line.graph(g)
g3
graph.isomorphic(g2, g3)
