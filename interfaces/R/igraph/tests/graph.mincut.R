
library(igraph)

g2 <- graph( c(1,2,2,3,3,4, 1,6,6,5,5,4, 4,1) )
E(g2)$capacity <- c(3,1,2, 10,1,3, 2)
graph.mincut(g2, value.only=FALSE)

