
library(igraph)

I <- matrix(sample(0:1, 35, replace=TRUE, prob=c(3,1)), nc=5)
g <- graph.incidence(I)
is.bipartite(g)$res

set.seed(42)
I <- matrix(sample(0:1, 35, replace=TRUE, prob=c(3,1)), nc=5)
g <- graph.incidence(I)
is.bipartite(g)

