
library(igraph)

g <- graph.formula( A-B-C, B-D )
are.connected(g, "A", "B")
are.connected(g, "B", "A")
are.connected(g, "A", "D")

g2 <- graph( c(1,2, 2,3, 3,4), dir=FALSE )
are.connected(g2, 1,2)
are.connected(g2, 3,2)
are.connected(g2, 4,1)

g3 <- graph.formula( A-+B-+C, B-+D )
are.connected(g3, "A", "C")
are.connected(g3, "A", "B")
are.connected(g3, "B", "A")

