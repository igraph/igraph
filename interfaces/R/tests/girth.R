
library(igraph)

## No circle in a tree
g <- graph.tree(1000, 3)
girth(g)

## The worst case running time is for a ring
g <- graph.ring(100)
girth(g)

