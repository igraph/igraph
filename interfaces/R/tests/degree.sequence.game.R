
library(igraph)

gc <- function(graph) {
  clu <- clusters(graph)
  induced.subgraph(graph, which(clu$membership==which.max(clu$csize)))
}

g <- gc(erdos.renyi.game(1000, 2/1000))

nG <- degree.sequence.game(degree(g), method="simple")
all(degree(nG) == degree(g))

nG <- degree.sequence.game(degree(g), method="vl")
all(degree(nG) == degree(g))
is.connected(nG)
is.simple(nG)

#####

g <- erdos.renyi.game(1000, 1/1000)

nG <- degree.sequence.game(degree(g), method="simple")
all(degree(nG) == degree(g))

g2 <- erdos.renyi.game(1000, 2/1000, dir=TRUE)

nG2 <- degree.sequence.game(degree(g, mode="out"),
                            degree(g, mode="in"),
                            method="simple")
all(degree(nG, mode="out") == degree(g, mode="out"))
all(degree(nG, mode="in") == degree(g, mode="in"))

