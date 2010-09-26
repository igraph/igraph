
library(igraph)

g <- graph.tree(7, 2)
g <- add.edges(g, c(2,7, 1,4))
g2 <- unfold.tree(g, roots=1)
g2
