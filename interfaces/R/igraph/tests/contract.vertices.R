
library(igraph)

set.seed(42)

g <- graph.ring(10)
g$name <- "Ring"
V(g)$name <- letters[1:vcount(g)]
E(g)$weight <- runif(ecount(g))
E(g)$weight

g2 <- contract.vertices(g, rep(1:5, each=2),
                        vertex.attr.comb=toString)

## graph and edge attributes are kept, vertex attributes are
## combined using the 'toString' function.
print(g2, g=TRUE, v=TRUE, e=TRUE)

