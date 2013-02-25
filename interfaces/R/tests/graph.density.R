
library(igraph)

g <- erdos.renyi.game(50, 4/50)
gd <- graph.density(g)
gd2 <- ecount(g) / vcount(g) / (vcount(g)-1) * 2
abs(gd - gd2) < 1e-14

####

g <- erdos.renyi.game(50, 4/50, dir=TRUE)
gd <- graph.density(g)
gd2 <- ecount(g) / vcount(g) / (vcount(g)-1)
abs(gd - gd2) < 1e-14

