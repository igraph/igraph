
library(igraph)

apl <- function(graph) {
  sp <- shortest.paths(graph, mode="out")
  if (is.directed(graph)) {
    diag(sp) <- NA
  } else {
    sp[lower.tri(sp, diag=TRUE)] <- NA
  }
  sp[sp=="Inf"] <- NA
  mean(sp, na.rm=TRUE)
}

giant.component <- function(graph, mode="weak") {
  clu <- clusters(graph, mode=mode)
  induced.subgraph(graph, which(clu$membership==which.max(clu$csize)))
}
  
g <- giant.component(erdos.renyi.game(100, 3/100))
abs(apl(g) - average.path.length(g)) < 1e-14

g <- giant.component(erdos.renyi.game(100, 6/100, dir=TRUE), mode="strong")
abs(apl(g) - average.path.length(g)) < 1e-14

g <- erdos.renyi.game(100, 2/100)
abs(apl(g) - average.path.length(g)) < 1e-14

g <- erdos.renyi.game(100, 4/100, dir=TRUE)
abs(apl(g) - average.path.length(g)) < 1e-14
