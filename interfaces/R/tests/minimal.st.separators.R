
library(igraph)

g <- graph.famous("Zachary")
msts <- minimal.st.separators(g)
all(sapply(msts, is.separator, graph=g))
