
library(igraph)

g <- erdos.renyi.game(50, 0.8)
livs <- largest.independent.vertex.sets(g)
unique(sapply(livs, length)) == independence.number(g)

ec <- sapply(seq_along(livs), function(x)
             ecount(induced.subgraph(g, livs[[x]])))
unique(ec)
