
library(igraph) ; igraph.options(print.full=TRUE)

g <- graph.formula(A-B-C-D-E)
g2 <- add.vertices(g, 4)
g2

g3 <- add.vertices(g, 3, attr=list(name=c("F","G","H"), weight=1:3))
g3
V(g3)
V(g3)$weight

