
context("Graphlets")

test_that("Getting subcliques works", {
  library(igraph)
  set.seed(42*42)
  g <- erdos.renyi.game(10, 4/10)
  E(g)$weight <- as.double(sample(1:10, ecount(g), replace=TRUE))
  ids <- 1:vcount(g)

  cl <- maximal.cliques(g)
  cl <- lapply(cl, "-", 1)[c(9, 2, 3, 10, 5, 7, 6, 1, 4, 8)]

  res <- .Call("R_igraph_subclique_next", g, E(g)$weight, ids, cl,
               PACKAGE="igraph")

  for (i in seq_along(res$graphs)) {
    V(res$graphs[[i]])$name <- res$ids[[i]]
    E(res$graphs[[i]])$weight <- res$weights[[i]]    
  }
  
  expect_that(res$thr, equals(c(7,2,5,3,1,7,3,2,4,7)))
  expect_that(res$next_thr, equals(c(Inf, 4, 8, 8, Inf, 9, 5, 3, 5, Inf)))
  expect_that(res$weights, equals(list(numeric(), c(4,8), c(8,8), c(9,8),
                                       numeric(), c(9,9), c(7,7,9,10,5),
                                       c(7,3,4,5), c(5,7,5,9,10),
                                       numeric())))
  expect_that(res$ids, equals(list(integer(), c(5,9,10), c(1,9,10),
                                   c(1,8,9), integer(), c(3,7,6),
                                   c(3,6,5,4), c(3,5,10,2), c(2,5,3,4),
                                   integer())))
  expect_that(sapply(res$graphs, vcount), equals(sapply(res$ids, length)))
  expect_that(sapply(res$graphs, ecount),
              equals(sapply(res$weights, length)))
})

sortgl <- function(x) {
  cl <- lapply(x$cliques, sort)
  n <- sapply(cl, length)
  list(cliques=cl[order(n)], thresholds=x$thresholds[order(n)])
}

test_that("Graphlets work for some simple graphs", {
  library(igraph)

  g <- graph.full(5)
  E(g)$weight <- 1
  gl <- graphlets.candidate.basis(g)

  expect_that(names(gl), equals(c("cliques", "thresholds")))
  expect_that(length(gl$cliques), equals(1))
  expect_that(sort(gl$cliques[[1]]), equals(1:vcount(g)))
  expect_that(gl$thresholds, equals(1))

  g2 <- graph.full(5)
  E(g2)$weight <- 1
  E(g2)[1%--%2]$weight <- 2
  gl2 <- sortgl(graphlets.candidate.basis(g2))

  expect_that(gl2, equals(list(cliques=list(1:2, 1:5), thresholds=c(2,1))))
})

test_that("Graphlets filtering works", {
  library(igraph)
  gt <- data.frame(from  =c("A", "A", "B", "B", "B", "C", "C", "D"),
                   to    =c("B", "C", "C", "D", "E", "D", "E", "E"),
                   weight=c( 8 ,  8 ,  8 ,  5 ,  5 ,  5 ,  5 ,  5 ))

  g <- graph.data.frame(gt, directed=FALSE, vertices=data.frame(LETTERS[1:5]))
  gl <- sortgl(graphlets.candidate.basis(g))

  expect_that(gl$cliques, equals(list(1:3, 2:5)))
  expect_that(gl$thresholds, equals(c(8, 5)))
})

## Naive version of graphlets

threshold.net <- function(graph, level) {
  N <- vcount(graph)
  graph.t <- delete.edges(graph, which(E(graph)$weight < level))

  clqt <- maximal.cliques(graph.t)
  clqt <- lapply(clqt, sort)
  clqt[order(sapply(clqt, length), decreasing=TRUE)]
}

graphlets.old <- function(graph) {

  if (!is.weighted(graph)) { stop("Graph not weighted") }
  if (min(E(graph)$weight) <= 0 || !is.finite(E(graph)$weight)) {
    stop("Edge weights must be non-negative and finite")
  }

  ## Do all thresholds
  cl <- lapply(sort(unique(E(graph)$weight)), function(w) {
    threshold.net(graph, w)
  })

  ## Put the cliques in one long list
  clv <- unlist(cl, recursive=FALSE)

  ## Sort the vertices within the cliques
  cls <- lapply(clv, sort)

  ## Delete duplicate cliques
  clu <- unique(cls)

  ## Delete cliques that consist of single vertices
  clf <- clu[sapply(clu, length) != 1]

  clf
}

test_that("Graphlets work for a bigger graph", {
  library(igraph)
  set.seed(42)
  g <- graph.famous("zachary")
  E(g)$weight <- sample(1:5, ecount(g), replace=TRUE)

  gl <- graphlets.candidate.basis(g)
  gl2 <- graphlets.old(g)

  glo <- sort(sapply(gl$cliques, paste, collapse="-"))
  gl2o <- sort(sapply(gl2, paste, collapse="-"))

  expect_that(glo, equals(gl2o))
})

graphlets.project.old <- function(graph, cliques, iter, Mu=NULL) {

  if (!is.weighted(graph)) { stop("Graph not weighted") }
  if (min(E(graph)$weight) <= 0 || !is.finite(E(graph)$weight)) {
    stop("Edge weights must be non-negative and finite")
  }
  if (length(iter) != 1 || !is.numeric(iter) ||
      !is.finite(iter) || iter != as.integer(iter)) {
    stop("`iter' must be a non-negative finite integer scalar")
  }

  clf <- cliques

  ## Create vertex-clique list first
  vcl <- vector(length=vcount(graph), mode="list")
  for (i in 1:length(clf)) {
    for (j in clf[[i]]) {
      vcl[[j]] <- c(vcl[[j]], i)
    }
  }

  ## Create edge-clique list from this, it is useful to have the edge list
  ## of the graph at hand
  el <- get.edgelist(graph, names=FALSE)
  ecl <- vector(length=ecount(graph), mode="list")
  for (i in 1:ecount(graph)) {
    edge <- el[i,]
    ecl[[i]] <- intersect(vcl[[edge[1]]], vcl[[edge[2]]])
  }

  ## We will also need a clique-edge list, the edges in the cliques
  system.time({
    cel <- vector(length=length(clf), mode="list")
    for (i in 1:length(ecl)) {
      for (j in ecl[[i]]) {
        cel[[j]] <- c(cel[[j]], i)
      }
    }
  })

  ## OK, we are ready to do the projection now
  if (is.null(Mu)) { Mu <- rep(1, length(clf)) }
  origw <- E(graph)$weight
  w <- numeric(length(ecl))
  a <- sapply(clf, function(x) length(x) * (length(x) + 1) / 2)
  for (i in 1:iter) {
    for (j in 1:length(ecl)) {
      w[j] <- sum(Mu[ ecl[[j]] ])
    }
    for (j in 1:length(clf)) {
      Mu[j] <- Mu[j] * sum(origw[cel[[j]]] / (w[cel[[j]]] + .0001)) / a[j]
    }
  }

  ## Sort the cliques according to their weights
  Smb <- sort(Mu, decreasing=TRUE, index=TRUE)
  list(cliques=clf[Smb$ix], Mu=Mu[Smb$ix])
}

test_that("Graphlet projection works", {
  library(igraph)

  D1 <- matrix(0, 5, 5)
  D2 <- matrix(0, 5, 5)
  D3 <- matrix(0, 5, 5)
  D1[1:3, 1:3] <- 2
  D2[3:5, 3:5] <- 3
  D3[2:5, 2:5] <- 1
  g <- graph.adjacency(D1 + D2 + D3, mode="undirected", weighted=TRUE)
  g <- simplify(g)

  gl <- graphlets.candidate.basis(g)
  glp <- graphlets(g)
  glp2 <- graphlets.project.old(g, cliques=gl$cliques, iter=1000)

  expect_that(glp, equals(glp2))
})
