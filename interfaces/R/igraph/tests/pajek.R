
library(igraph)

g <- graph.ring(9)
V(g)$color <- c("red", "green", "yellow")

tc <- rawConnection(raw(0), "w")
write.graph(g, format="pajek", file=tc)
cat(rawToChar(rawConnectionValue(tc)))
close(tc)

