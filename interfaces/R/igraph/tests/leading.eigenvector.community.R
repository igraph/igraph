
library(igraph)

## Check-test

f <- function(membership, community, value, vector, multiplier, extra) { 
  M <- sapply(1:length(vector), function(x) {
    v <- rep(0, length(vector))
    v[x] <- 1
    multiplier(v)
  })
  ev <- eigen(M)
  ret <- 0
  if (abs(ev$values[1] - value) > 1e-10) { ret <- 1 }
  if (sign(ev$vectors[1,1]) != sign(vector[1])) { ev$vectors <- -ev$vectors }
  if (any(abs(ev$vectors[,1] - vector) > 1e-10)) { ret <- 1 }

  if (ret) {
    print("Error!")
    print(M)
    print(c(value, ev$values[1]))
    print(cbind(vector, ev$vectors[,1]))
  }
  
  ret
}

g <- graph.famous("Zachary")
lc <- leading.eigenvector.community(g, callback=f)
lc
lc$modularity == modularity(g, lc$membership)
membership(lc)
length(lc)
sizes(lc)

## Stress-test

for (i in 1:100) {
  g <- erdos.renyi.game(20, sample(5:40, 1), type="gnm")
  lec1 <- leading.eigenvector.community(g)
  lec2 <- leading.eigenvector.community(g)
  if (length(membership(lec1)) != length(membership(lec2)) ||
      any(membership(lec1) != membership(lec2))) {
    print("Foobar!")
    break
  }
}

