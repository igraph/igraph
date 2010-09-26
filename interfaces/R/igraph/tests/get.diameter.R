
library(igraph)

g <- graph.ring(10)
E(g)$weight <- sample(seq_len(ecount(g)))
d <- diameter(g)
gd <- get.diameter(g)
sp <- shortest.paths(g)

d == max(sp)
sp[ gd[1], gd[length(gd)] ] == d

d <- diameter(g, weights=NA)
gd <- get.diameter(g, weights=NA)
sp <- shortest.paths(g, weights=NA)

d == max(sp)
length(gd) == d + 1
sp[ gd[1], gd[length(gd)] ] == d

