
library(igraph)

g <- graph.famous("Zachary")
oc <- optimal.community(g)
oc
membership(oc)
abs(modularity(g, oc$membership) - oc$modularity) < 1e-14
length(oc)
sizes(oc)
