
library(igraph) ; igraph.options(print.full=TRUE)

set.seed(42)

g <- graph.full.bipartite(10,5)
proj <- bipartite.projection(g)
graph.isomorphic(proj[[1]], graph.full(10))
graph.isomorphic(proj[[2]], graph.full(5))

M <- matrix(0, nr=5, nc=3)
rownames(M) <- c("Alice", "Bob", "Cecil", "Dan", "Ethel")
colnames(M) <- c("Party", "Skiing", "Badminton")
M[] <- sample(0:1, length(M), replace=TRUE)
M
g2 <- graph.incidence(M)
g2$name <- "Event network"
proj2 <- bipartite.projection(g2)
print(proj2[[1]], g=TRUE, e=TRUE)
print(proj2[[2]], g=TRUE, e=TRUE)

bipartite.projection.size(g2)
