
library(igraph)

E <- rbind( c(1,3,3), c(3,4,1), c(4,2,2), c(1,5,1), c(5,6,2), c(6,2,10))
colnames(E) <- c("from", "to", "capacity")
g1 <- graph.data.frame(as.data.frame(E))
graph.maxflow(g1, source="1", target="2")

