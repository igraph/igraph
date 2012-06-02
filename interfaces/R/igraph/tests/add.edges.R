
library(igraph) ; igraph.options(print.full=TRUE)

g <- graph.empty(10)
g2 <- add.edges(g, c(1,2, 2,3, 3,4, 1,6, 1,7, 9,10) )
g2

g3 <- add.edges(g, c(1,5, 2,6, 3,10, 4,5), attr=list(weight=c(1,2,1,-1)) )
g3
E(g3)$weight

g4 <- add.edges(g2, c(1,4, 4,6, 7,1), attr=list(weight=c(-1,1,-2.5)))
g4
E(g4)$weight

g5 <- add.edges(g3, c(10,9, 10,10, 1,1), attr=list(weight=c(100,100,100)) )
g5
E(g5)$weight


