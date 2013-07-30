

library(igraph)
set.seed(42)

for (i in 1:100) {
  cat(".")
  g <- erdos.renyi.game(100, 10/100)
  
  s1 <- scan1(g)
  
  s1a <- scan1.approx(g, 20)$res
  cor(s1, s1a)

  E <- graph.eigen(g, which=list(howmany=20, pos="LM"))
  s1aa <- scan1.approx.eigen(g, E$values, E$vectors)
  cor(s1, s1aa)

  E2 <- eigen(get.adjacency(g, sparse=FALSE))
  s1aaa <- colSums(E2$values ^3 * t(E2$vectors)^2 / 2) + degree(g)
  max(abs(s1aaa - s1))

  std <- function(x) {
    x <- zapsmall(x)
    apply(x, 2, function(col) {
      if (any(col < 0) && col[which(col != 0)[1]] < 0) { -col } else { col }
    })
  }

  ord <- order(abs(E2$values), decreasing=TRUE)
  max(abs(E$values - E2$values[ord[1:20]]))
  max(abs(std(E$vectors) - std(E2$vectors[,ord[1:20]])))

}
