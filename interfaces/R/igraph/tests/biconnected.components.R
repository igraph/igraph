
library(igraph)

g <- graph.full(5) %du% graph.full(5)
clu <- clusters(g)$membership
g <- add.edges(g, c(match(1,clu), match(2,clu)) )

biconnected.components(g)
