
library(igraph)

g <- graph.formula(a -- 1:3 -- 5 -- 2:4 -- b, 1 -- 2, 3 -- 4)
stsep <- minimal.st.separators(g)
sapply(stsep, is.minimal.separator, graph=g)

