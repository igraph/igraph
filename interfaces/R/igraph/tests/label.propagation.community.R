
library(igraph)

g <- graph.famous("Zachary")
set.seed(42)
lpc <- label.propagation.community(g)
lpc
lpc$modularity == modularity(g, lpc$membership)
membership(lpc)

