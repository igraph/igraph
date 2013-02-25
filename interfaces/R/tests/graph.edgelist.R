
library(igraph)

g <- erdos.renyi.game(50, 5/50)
el <- get.edgelist(g)
g2 <- graph.edgelist(el, dir=FALSE)
graph.isomorphic(g, g2)

####

g <- erdos.renyi.game(50, 5/50, dir=TRUE)
el <- get.edgelist(g)
g2 <- graph.edgelist(el, dir=TRUE)
graph.isomorphic(g, g2)

####

g <- erdos.renyi.game(26, 5/26, dir=TRUE)
el <- get.edgelist(g)
n <- letters[1:26]
names(n) <- 1:26
mode(el) <- "character"
el[] <- n[el]
g2 <- graph.edgelist(el, dir=TRUE)
graph.isomorphic(g, g2)

