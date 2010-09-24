
library(igraph)

set.seed(42)

## Some trivial ones
g <- graph.ring(10)
graph.knn(g)
g2 <- graph.star(10)
graph.knn(g2)

## A scale-free one, try to plot 'knnk'
g3 <- simplify(ba.game(1000, m=5))
graph.knn(g3)

## A random graph
g4 <- random.graph.game(1000, p=5/1000)
graph.knn(g4)

## A weighted graph
g5 <- graph.star(10)
E(g5)$weight <- seq(ecount(g5))
graph.knn(g5)
