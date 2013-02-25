
library(igraph)

g <- erdos.renyi.game(50,20/50)
lc <- largest.cliques(g)
unique(sapply(lc, function(x) graph.density(induced.subgraph(g, x))))
