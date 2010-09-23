
library(igraph)

edges <- matrix(c("s", "a", 2,
                  "s", "b", 4,
                  "a", "t", 4,
                  "b", "t", 2,
                  "a", "1", 1,
                  "a", "2", 1,
                  "a", "3", 2,
                  "1", "b", 1,
                  "2", "b", 2,
                  "3", "b", 1),
                byrow=TRUE, ncol=3,
                dimnames=list(NULL, c("from", "to", "weight")))
edges <- as.data.frame(edges)
edges[[3]] <- as.numeric(as.character(edges[[3]]))

g <- graph.data.frame(as.data.frame(edges))

get.all.shortest.paths(g, "s", "t", weights=NA)
get.all.shortest.paths(g, "s", "t")

## TODO

## E(g)$weight <- E(g)$weight - 1
## get.all.shortest.paths(g, "s", "t")

