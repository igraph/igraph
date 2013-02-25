
library(igraph)

g <- erdos.renyi.game(50, 0.8)
ivs <- independent.vertex.sets(g, min=independence.number(g))
ec <- sapply(seq_along(ivs), function(x)
             ecount(induced.subgraph(g, ivs[[x]])))
unique(ec)
