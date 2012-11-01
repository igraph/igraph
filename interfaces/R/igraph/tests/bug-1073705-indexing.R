
library(igraph)

g <- graph.ring(10)
g[1,2,attr="weight"] <- 0
g
g[]

el <- get.edgelist(g)
g[from=el[,1], to=el[,2], attr="sim"] <- rep(0:1, length=ecount(g))
g
g[attr="sim"]

V(g)$name <- letters[seq_len(vcount(g))]
el <- get.edgelist(g)
g[from=el[,1], to=el[,2], attr="sim"] <- rep(1:0, length=ecount(g))
g
g[attr="sim"]

