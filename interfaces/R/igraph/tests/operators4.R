
library(igraph) 

g1 <- g2 <- graph.ring(10)
g1$foo <- "bar"
V(g1)$name <- letters[ 1:10]
V(g2)$name <- letters[11:20]
E(g1)$weight <- 1:10
E(g2)$weight <- 10:1

V(g1)$a1 <- 1:10
V(g2)$a2 <- 11:20

E(g1)$b1 <- 1:10
E(g2)$b2 <- 11:20

g1 + g2

V(g1+g2)$name
V(g1+g2)$a1
V(g1+g2)$a2

E(g1+g2)$weight
E(g1+g2)$b1
E(g1+g2)$b2

