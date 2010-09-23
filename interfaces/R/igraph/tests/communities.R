
library(igraph)

set.seed(42)

F <- list("edge.betweenness.community", "fastgreedy.community",
          "label.propagation.community", "leading.eigenvector.community",
          "multilevel.community", "optimal.community", "spinglass.community",
          "walktrap.community")

karate <- graph.famous("Zachary")

for (f in F) {
  print(f)
  f <- get(f)
  comm <- f(karate)
  
  c1 <- abs(modularity(comm) - modularity(karate, membership(comm))) < 1e-9
  
  cc <- communities(comm)
  c2 <- all(!duplicated(unlist(cc)))
  
  c3 <- all(unlist(cc) <= vcount(karate) & unlist(cc) >= 1)

  c4 <- length(comm) == max(membership(comm))
  
  print(c(c1, c2, c3, c4))
}

fc <- fastgreedy.community(karate)
modularity(karate, cutat(fc, no=1))
modularity(karate, cutat(fc, no=2))
modularity(karate, cutat(fc, no=3))
modularity(karate, cutat(fc, no=4))

crossing(fc, karate)

