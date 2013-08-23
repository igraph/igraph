
context("Maximal cliques")

PP <- list()
XX <- list()

bk3 <- function(graph, min=0, max=Inf) {

  Gamma <- function(v) { neighbors(graph, v) }

  bkpivot <- function(P, R, X) {
    PP <<- c(PP, list(p1=sort(P)))
    XX <<- c(XX, list(x1=sort(X)))
    if (length(P) == 0 && length(X) == 0) {
      if (length(R) >= min && length(R) <= max) { list(R) } else { list() }
    } else if (length(P) != 0) { 
      psize <- sapply(c(P, X),
                      function(u) length(intersect(P, Gamma(u))))
      u <- c(P, X)[which.max(psize)]
      
      pres <- list()
      for (v in setdiff(P, Gamma(u))) {
        pres <- c(pres, bkpivot(intersect(P, Gamma(v)), c(R, v),
                                intersect(X, Gamma(v))))
        P <- P[P != v]
        X <- c(X, v)
        PP <<- c(PP, list(p2=sort(P)))
        XX <<- c(XX, list(x2=sort(X)))
      }
      pres
    }
  }

  res <- list()
  cord <- order(graph.coreness(graph))
  for (v in seq_along(cord)) {
    if (v != length(cord)) {
      P <- intersect(Gamma(cord[v]), cord[(v+1):length(cord)])
    } else {
      P <- numeric()
    }
    if (v != 1) {
      X <- intersect(Gamma(cord[v]), cord[1:(v-1)])
    } else {
      X <- numeric()
    }
    res <- c(res, bkpivot(P, cord[v], X))
  }
  res
}

mysort <- function(x) {
  xl <- sapply(x, length)
  x <- lapply(x, sort)
  xc <- sapply(x, paste, collapse="-")
  x[order(xl, xc)]
}

###################################################################

test_that("Maximal cliques work", {
  library(igraph)
  set.seed(42)
  G <- erdos.renyi.game(1000, 1000, type="gnm")
  cli <- graph.full(10)
  for (i in 1:10) {
    G <- permute.vertices(G, sample(vcount(G)))
    G <- G %u% cli
  }
  G <- simplify(G)

  cl1 <- mysort(bk3(G, min=3))
  cl2 <- mysort(maximal.cliques(G, min=3))

  expect_that(cl1, is_identical_to(cl2))
})

###################################################################

PP2 <- list()
XX2 <- list()

bk4 <- function(graph, min=0, max=Inf) {

  Gamma <- function(v) { neighbors(graph, v) }

  bkpivot <- function(PX, R) {
    P <- if (PX$PE >= PX$PS) { PX$PX[PX$PS:PX$PE] } else { numeric() }
    X <- if (PX$XE >= PX$XS) { PX$PX[PX$XS:PX$XE] } else { numeric() }
    PP2 <<- c(PP2, list(p1=sort(P)))
    XX2 <<- c(XX2, list(x1=sort(X)))
    if (length(P) == 0 && length(X) == 0) {
      if (length(R) >= min && length(R) <= max) { list(R) } else { list() }
    } else if (length(P) != 0) {
      psize <- sapply(c(P, X), function(u)
                      length(intersect(P, Gamma(u))))
      u <- c(P, X)[which.max(psize)]

      pres <- list()
      for (v in setdiff(P, Gamma(u))) {
        
        p0 <- if (PX$PS > 1) { PX$PX[1:(PX$PS-1)] } else { numeric() }
        p1 <- setdiff(P, Gamma(v))
        p2 <- intersect(P, Gamma(v))
        x1 <- intersect(X, Gamma(v))
        x2 <- setdiff(X, Gamma(v))
        x0 <- if (PX$XE < length(PX$PX)) {
          PX$PX[(PX$XE+1):length(PX$PX)]
        } else {
          numeric()
        }

        newPX <- list(PX=c(p0, p1, p2, x1, x2, x0),
                      PS=length(p0) + length(p1) + 1,
                      PE=length(p0) + length(p1) + length(p2),
                      XS=length(p0) + length(p1) + length(p2) + 1,
                      XE=length(p0) + length(p1) + length(p2) + length(x1))
        
        pres <- c(pres, bkpivot(newPX, c(R, v)))

        vpos <- which(PX$PX==v)
        tmp <- PX$PX[PX$PE]
        PX$PX[PX$PE] <- v
        PX$PX[vpos] <- tmp
        PX$PE <- PX$PE - 1
        PX$XS <- PX$XS - 1
        P <- if (PX$PE >= PX$PS) { PX$PX[PX$PS:PX$PE] } else { numeric() }
        X <- if (PX$XE >= PX$XS) { PX$PX[PX$XS:PX$XE] } else { numeric() }
        PP2 <<- c(PP2, list(p2=sort(P)))
        XX2 <<- c(XX2, list(x2=sort(X)))
        if (any(duplicated(PX$PX))) { stop("foo2") }
      }
      pres
    }
  }

  res <- list()
  cord <- order(graph.coreness(graph))
  for (v in seq_along(cord)) {
    if (v != length(cord)) {
      P <- intersect(Gamma(cord[v]), cord[(v+1):length(cord)])
    } else {
      P <- numeric()
    }
    if (v != 1) {
      X <- intersect(Gamma(cord[v]), cord[1:(v-1)])
    } else {
      X <- numeric()
    }
    PX <- list(PX=c(P, X), PS=1, PE=length(P),
               XS=length(P)+1, XE=length(P)+length(X))
    res <- c(res, bkpivot(PX, cord[v]))
  }
  res    
}

#################################################################

test_that("Maximal cliques work", {
  library(igraph)
  set.seed(42)
  G <- erdos.renyi.game(1000, 1000, type="gnm")
  cli <- graph.full(10)
  for (i in 1:10) {
    G <- permute.vertices(G, sample(vcount(G)))
    G <- G %u% cli
  }
  G <- simplify(G)

  cl1 <- mysort(bk4(G, min=3))
  cl2 <- mysort(maximal.cliques(G, min=3))

  expect_that(cl1, is_identical_to(cl2))
})
