
library(igraph)

g <- erdos.renyi.game(100, 2/100)
g2 <- as.directed(g, mode="mutual")
g3 <- as.directed(g, mode="arbitrary")

all(degree(g) == degree(g3))
all(degree(g) == degree(g2) / 2)

graph.isomorphic(g, as.undirected(g2))
graph.isomorphic(g, as.undirected(g3))

ge <- graph.empty(10)
summary(as.directed(ge))
summary(as.directed(ge, "arbitrary"))

ge2 <- graph.empty()
summary(as.directed(ge2))
summary(as.directed(ge2, "arbitrary"))

#### Keeping attributes

g <- graph.formula( A-B-C, D-A, E )
g$name <- "Small graph"
g2 <- as.directed(g, mode="mutual")
g3 <- as.directed(g, mode="arbitrary")
g2
g3
g2$name
g3$name

set.seed(42)
E(g)$weight <- runif(ecount(g))
print(as.directed(g, "mutual"), e=TRUE)
print(as.directed(g, "arbitrary"), e=TRUE)

