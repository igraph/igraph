
## Algorithm 1:
## D is NxN non-negative weighted matrix
## level is a number between minimum and maximum element in D
## puts a Threshold on D at level and returns all maximal cliques of
## the binary network 

threshold.net <- function(graph, level) {
  N <- vcount(graph)
  graph.t <- delete.edges(graph, which(E(graph)$weight < level))

  clqt <- maximal.cliques(graph.t)
  clqt <- lapply(clqt, sort)
  clqt[order(sapply(clqt, length), decreasing=TRUE)]
}


## Algorithm 2:
## D is NxN non negative weighted matrix
##  projects the network on the basis elements(Bc) and finds Mu

project.net <- function(graph, Bc, iter) {
  K <- length(Bc)
  a <- sapply(Bc, function(x) length(x) * (length(x)+1) / 2)
  Mu <- rep(1, length(Bc))

  sele <- lapply(1:K, function(j) {
    unique(as.vector(graph[ Bc[[j]], Bc[[j]], edges=TRUE, sparse=FALSE ]))
  })
  origw <- lapply(1:K, function(j) {
    E(graph)$weight[sele[[j]]]
  })
  for (i in 1:iter) {
    w <- numeric(ecount(graph))
    for (j in 1:K) { w[sele[[j]]] <- w[sele[[j]]] + Mu[j] }
    for (j in 1:K) {
      Qt <- Bc[[j]]
      Mu[j] <- Mu[j] * sum(origw[[j]] / (w[sele[[j]]] + .0001)) / a[j]
    }
  }
  Mu
}


## Main:
## Finds all basis(Bc) for the network D (using algorithm 1), projects
## on the Bc and finds Muc (using algorithm 2). Input D must be
## symettric with all nonnegative elements

graphlets <- function(graph, iter) {  

  if (!is.weighted(graph)) { stop("Graph not weighted") }
  if (min(E(graph)$weight) <= 0 || !is.finite(E(graph)$weight)) {
    stop("Edge weights must be non-negative and finite")
  }
  if (length(iter) != 1 || !is.numeric(iter) ||
      !is.finite(iter) || iter != as.integer(iter)) {
    stop("`iter' must be a non-negative finite integer scalar")
  }

  ## This is slow, but readable, will speed it up later
  L <- c(unique(E(graph)$weight), 0.0001)
  Bc <- list()
  for (level in L) {
    Bt <- threshold.net(graph, level)
    Bt <- Bt[sapply(Bt, length) > 1]
    Bc <- unique(c(Bc, Bt))
  }

  Mu <- project.net(graph, Bc, iter)
  Smb <- sort(Mu, decreasing=TRUE, index=TRUE)

  list(Bc=Bc[Smb$ix], Muc=Mu[Smb$ix])
}

#################
## Example code

function() {
  library(igraph)
  D1 <- matrix(0, 5, 5)
  D2 <- matrix(0, 5, 5)
  D3 <- matrix(0, 5, 5)
  D1[1:3, 1:3] <- 2
  D2[3:5, 3:5] <- 3
  D3[2:5, 2:5] <- 1
  
  g <- graph.adjacency(D1 + D2 + D3, mode="undirected", weighted=TRUE)
  
  gl <- graphlets(g, iter=1000)

  g <- simplify(g)
  V(g)$color <- "white"
  E(g)$label <- E(g)$weight
  E(g)$label.cex <- 2
  E(g)$color <- "black"

  plot.new()
  layout(matrix(1:6, nrow=2, byrow=TRUE))
  co <- layout.kamada.kawai(g)
  par(mar=c(1,1,1,1))
  plot(g, layout=co)
  for (i in 1:length(gl$Bc)) {
    sel <- gl$Bc[[i]]
    V(g)$color <- "white"
    V(g)[sel]$color <- "#E495A5"
    E(g)$width <- 1
    E(g)[ V(g)[sel] %--% V(g)[sel] ]$width <- 2
    E(g)$label <- ""
    E(g)[ width == 2 ]$label <- round(gl$Muc[i], 2)
    E(g)$color <- "black"
    E(g)[ width == 2 ]$color <- "#E495A5"
    plot(g, layout=co)
  }
}

