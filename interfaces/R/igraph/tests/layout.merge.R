
library(igraph)

set.seed(42)

g <- list(graph.ring(10), graph.ring(5))
l <- lapply(g, layout.mds)
l

layout.merge(g, l)

##########
# Stress test

for (i in 1:100) {
  g <- erdos.renyi.game(100, 2/100)
  l <- layout.mds(g)
}

