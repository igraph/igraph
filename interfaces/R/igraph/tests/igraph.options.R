
library(igraph)

g <- graph.ring(5)
V(g)$name <- letters[1:5]
E(g)$weight <- 5:1
g$name <- "A small ring"

g

igraph.options(print.vertex.attributes=TRUE)
getIgraphOpt("print.vertex.attributes")
g

igraph.options(print.edge.attributes=TRUE, print.vertex.attributes=FALSE)
g

igraph.options(print.graph.attributes=TRUE)
g
