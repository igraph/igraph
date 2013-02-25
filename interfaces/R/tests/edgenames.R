
library(igraph) ; igraph.options(print.full=TRUE)

## named edges
igraph.options(print.edge.attributes = TRUE)
g <- graph.ring(10)
E(g)$name <- letters[1:ecount(g)]
delete.edges(g, c("b", "d", "e"))

## named vertices
g <- graph.ring(10)
V(g)$name <- letters[1:vcount(g)]
delete.edges(g, c("a|b", "f|g", "c|b"))

## no names at all, but select edges based on vertices
g <- graph.ring(10)
delete.edges(g, c("1|2", "8|7", "1|10"))

## mix edge names and vertex names
g <- graph.ring(10)
V(g)$name <- letters[1:vcount(g)]
E(g)$name <- LETTERS[1:ecount(g)]
delete.edges(g, c("a|b", "F", "j|i"))
