
library(igraph)

set.seed(42)

gc <- function(graph) {
  cl <- clusters(graph)
  induced.subgraph(graph, which(cl$membership==which.max(cl$csize)))
}

rg <- function(n) {
  gc(erdos.renyi.game(n, 1/n))
}
  
G <- lapply(1:30, function(x) rg(sample(100, 1)))
Gsize <- sapply(G, vcount)

allg <- graph.disjoint.union(G)
clu <- clusters(allg)
all(table(clu$membership) == clu$csize)
all(sort(clu$csize)==sort(Gsize))
clu$no == length(G)

