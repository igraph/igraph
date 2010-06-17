
library(igraph)

g <- graph.famous("Zachary")
fc <- fastgreedy.community(g)
fc
modularity(g, fc$membership) == max(fc$modularity)
membership(fc)
length(fc)
sizes(fc)
