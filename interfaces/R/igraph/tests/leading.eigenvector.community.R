
library(igraph)

g <- graph.famous("Zachary")
lc <- leading.eigenvector.community(g)
lc
lc$modularity == modularity(g, lc$membership)
membership(lc)
length(lc)
sizes(lc)

## Stress-test

for (i in 1:1000) {
  g <- erdos.renyi.game(20, sample(5:40, 1), type="gnm")
  lec1 <- leading.eigenvector.community(g)
  lec2 <- leading.eigenvector.community(g)
  if (length(membership(lec1)) != length(membership(lec2)) ||
      any(membership(lec1) != membership(lec2))) {
    print("Foobar!")
    break
  }
}

