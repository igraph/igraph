
library(igraph)

data <- read.table("bug-1019624.txt")
data

g <- graph.adjacency(as.matrix(data), weighted=TRUE)
g
E(g)$weight
str(g)

