
library(igraph)

adj <- matrix(1, nrow=11, ncol=11) - diag(11)
g <- graph.adjacency(adj)
suppressWarnings(largest.cliques(g))

