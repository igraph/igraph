
library(igraph)
options(width=60)

kar <- graph.famous("zachary")

fc <- fastgreedy.community(kar)

modularity(kar, membership(fc))
modularity(kar, membership(fc), weights=rep(1, ecount(kar)))

B1 <- mod.matrix(kar, membership(fc))
B2 <- mod.matrix(kar, membership(fc), weights=rep(1, ecount(kar)))

max(abs(B1 - B2))

