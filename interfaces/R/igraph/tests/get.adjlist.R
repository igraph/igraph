
library(igraph)

g <- erdos.renyi.game(50, 2/50)
al <- get.adjlist(g)
g2 <- graph.adjlist(al, mode="all")
graph.isomorphic(g, g2)
graph.isomorphic.vf2(g, g2, vertex.color1=1:vcount(g),
                     vertex.color2=1:vcount(g2))

####

el <- get.adjedgelist(g)
for (i in 1:vcount(g)) {
  a <- as.numeric(E(g)[adj(i)])
  if (length(a) != length(el[[i]]) ||
      any(sort(el[[i]]) != sort(a))) { print("Foobar!"); break; }
}

g <- erdos.renyi.game(50, 4/50, directed=TRUE)
el1 <- get.adjedgelist(g, mode="out")
el2 <- get.adjedgelist(g, mode="in")
for (i in 1:vcount(g)) {
  a <- as.numeric(E(g)[from(i)])
  if (length(a) != length(el1[[i]]) ||
      any(sort(el1[[i]]) != sort(a))) { print("Foobar!"); break; }
}
for (i in 1:vcount(g)) {
  a <- as.numeric(E(g)[to(i)])
  if (length(a) != length(el2[[i]]) ||
      any(sort(el2[[i]]) != sort(a))) { print("Foobar!"); break; }
}
