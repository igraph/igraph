
library(igraph) ; igraph.options(print.full=TRUE)

g <- graph.formula(A:B:C - D:E:F, D-E-F)
g2 <- delete.edges(g, E(g, P=c("D", "E")))
g2
