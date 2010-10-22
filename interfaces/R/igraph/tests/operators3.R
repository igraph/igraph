
library(igraph)

g <- graph.ring(10)
V(g)$name <- letters[1:10]
E(g)$name <- LETTERS[1:10]

g <- g - c("a", "b")
g

g <- g - edge("e|f")
g

g <- g - edge("H")
g

g <- graph.ring(10)
V(g)$name <- letters[1:10]

g <- g - path("a", "b", "c", "d")
g

g - V(g)[c('d', 'g')]

g - E(g)['f' %--% 'g']
