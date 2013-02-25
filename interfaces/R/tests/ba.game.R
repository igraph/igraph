
# sink("ba.game.Rout.save")

library(igraph) ; igraph.options(print.full=TRUE)

g <- ba.game(100, m=2)
ecount(g)
vcount(g)
is.simple(g)

g2 <- ba.game(100, m=2, algorithm="psumtree-multiple")
ecount(g2)
vcount(g2)
is.simple(g2)

g3 <- ba.game(100, m=2, algorithm="bag")
ecount(g3)
vcount(g3)
is.simple(g3)

set.seed(1234)
g4 <- ba.game(10, m=1, algorithm="bag", start.graph=graph.empty(5))
ecount(g4)
vcount(g4)
degree(g4)

g6 <- ba.game(10, m=1, algorithm="bag", start.graph=graph.star(10))
g6

g7 <- ba.game(10, m=3, algorithm="psumtree-multiple",
              start.graph=graph.empty(5))
g7

g8 <- ba.game(10, m=3, algorithm="psumtree-multiple",
              start.graph=graph.star(5))
g8

g9 <- ba.game(10, m=3, algorithm="psumtree-multiple",
              start.graph=graph.star(10))
g9

g10 <- ba.game(10, m=3, start.graph=graph.empty(5))
g10

g11 <- ba.game(10, m=3, start.graph=graph.star(5))
g11

g12 <- ba.game(10, m=3, start.graph=graph.star(10))
g12

# sink(NULL)
