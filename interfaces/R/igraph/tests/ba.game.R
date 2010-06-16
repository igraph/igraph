
# sink("ba.game.Rout.save")

library(igraph)

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

# sink(NULL)
