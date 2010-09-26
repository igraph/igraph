
library(igraph)

set.seed(42)
g <- graph.formula(A+-+B, A--+C, C+-+D)
g$name <- "Tiny graph"
E(g)$weight <- rnorm(ecount(g))

g2 <- as.undirected(g, mode="collapse")
g3 <- as.undirected(g, mode="each")
g4 <- as.undirected(g, mode="mutual")

print(g, e=TRUE, g=TRUE)
print(g2, e=TRUE, g=TRUE)
print(g3, e=TRUE, g=TRUE)
print(g4, e=TRUE, g=TRUE)

