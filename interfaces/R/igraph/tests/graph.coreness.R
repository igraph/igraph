
library(igraph)

g <- graph.ring(10)
g <- add.edges(g, c(1,2, 2,3, 1,3))
graph.coreness(g)               

