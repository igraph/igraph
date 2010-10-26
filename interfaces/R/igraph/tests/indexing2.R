
library(igraph)

g <- graph.empty(10)
g[1,2] <- TRUE
g
g[2,1] <- TRUE
g
g[2,1] <- NULL
g
g[1,2] <- FALSE
g

g <- graph.empty(10)
g[-1,1] <- 1
g

## Weighted graph

igraph.options(print.edge.attributes=TRUE)
g <- graph.empty(10)
g <- set.edge.attribute(g, "weight", c(), 1)
g[1,2] <- 1
g
g[2,1] <- 2
g
g[1,2] <- 3
g
g[1:2,2:3] <- -1
g
g[1,2] <- NULL
g

## Using vertex names

g <- graph.empty(10)
V(g)$name <- letters[1:vcount(g)]
g['a','b'] <- TRUE
g['b','c'] <- TRUE
g
g[c('a','f'),c('f','a')] <- TRUE
g
g[c('a','c','h'), c('a','a','a'), attr="weight"] <- 3
g

