
library(igraph)

g <- erdos.renyi.game(1000, 1/1500)
G <- decompose.graph(g)
clu <- clusters(g)
Gsizes <- sapply(G, vcount)
all(sort(clu$csize) == sort(Gsizes))
