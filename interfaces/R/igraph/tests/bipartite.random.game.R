
library(igraph)

set.seed(42)
g1 <- bipartite.random.game(10, 5, type="gnp", p=.1)
str(g1)

g2 <- bipartite.random.game(10, 5, type="gnp", p=.1, directed=TRUE)
str(g2)

g3 <- bipartite.random.game(10, 5, type="gnp", p=.1, directed=TRUE, mode="in")
str(g3)

g4 <- bipartite.random.game(10, 5, type="gnm", m=8)
str(g4)

g5 <- bipartite.random.game(10, 5, type="gnm", m=8, directed=TRUE)
str(g5)

g6 <- bipartite.random.game(10, 5, type="gnm", m=8, directed=TRUE, mode="in")
str(g6)

#####

library(igraph)

g7 <- bipartite.random.game(10, 5, type="gnp", p=0.9999, directed=TRUE,
                            mode="all")
str(g7)

g8 <- bipartite.random.game(10, 5, type="gnm", m=99, directed=TRUE,
                            mode="all")
str(g8)


