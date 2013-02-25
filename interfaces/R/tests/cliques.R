
library(igraph)

set.seed(42)

check.clique <- function(graph, vids) {
  s <- induced.subgraph(graph, vids)
  ecount(s) == vcount(s) * (vcount(s)-1) / 2
}

g <- erdos.renyi.game(100, 0.3)
clique.number(g)
sapply(cliques(g, min=6), check.clique, graph=g)
sapply(largest.cliques(g), check.clique, graph=g)

## To have a bit less maximal cliques, about 100-200 usually
g <- erdos.renyi.game(100, 0.03)
all(sapply(maximal.cliques(g), check.clique, graph=g))

