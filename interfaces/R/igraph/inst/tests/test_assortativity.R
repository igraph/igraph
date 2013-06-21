
context("assortativity")

test_that("assortativity works", {
  library(igraph)
  
  g <- read.graph(f <- gzfile("celegansneural.gml.gz"), format="gml")

  assR <- function(graph) { 
    indeg <- degree(graph, mode="in")
    outdeg <- degree(graph, mode="out")
    el <- get.edgelist(graph, names=FALSE)
    J <- outdeg[el[,1]]-1
    K <- indeg[el[,2]]-1
    num <- sum(J*K) - sum(J)*sum(K)/ecount(graph)
    den1 <- sum(J*J) - sum(J)^2/ecount(graph)
    den2 <- sum(K*K) - sum(K)^2/ecount(graph)
    num / sqrt(den1) / sqrt(den2)
  }

  asd <- assortativity.degree(g)
  as <- assortativity(g, degree(g, mode="out"), degree(g, mode="in"))
  as2 <- assR(g)

  expect_that(asd, equals(as))
  expect_that(asd, equals(as2))

  asu <- assortativity.degree(simplify(as.undirected(g, mode="collapse")))
  expect_that(asu, equals(-0.16319921031570466807))

  p <- read.graph(f <- gzfile("power.gml.gz"), format="gml")
  p.asd <- assortativity.degree(p)
  p.as <- assortativity(p, degree(p))
  p.as2 <- assR(as.directed(p, mode="mutual"))

  expect_that(p.asd, equals(p.as))
  expect_that(p.asd, equals(p.as2))
})

test_that("nominal assortativity works", {
  library(igraph)
  
  o <- read.graph(f <- gzfile("football.gml.gz"), format="gml")
  o <- simplify(o)
  an <- assortativity.nominal(o, V(o)$value+1)

  el <- get.edgelist(o, names=FALSE)
  etm <- matrix(0, nr=max(V(o)$value)+1, nc=max(V(o)$value)+1)
  for (e in 1:nrow(el)) {
    t1 <- V(o)$value[ el[e,1] ]+1
    t2 <- V(o)$value[ el[e,2] ]+1
    etm[t1, t2] <- etm[t1, t2] + 1
    etm[t2, t1] <- etm[t2, t1] + 1
  }
  etm <- etm/sum(etm)
  an2 <- ( sum(diag(etm))-sum(etm %*% etm) ) / ( 1-sum(etm %*% etm) )

  expect_that(an, equals(an2))
})
