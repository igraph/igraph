

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

## Naive algorithm, do all thresholds for the whole graph

graphlets <- function(graph, iter) {

  if (!is.weighted(graph)) { stop("Graph not weighted") }
  if (min(E(graph)$weight) <= 0 || !is.finite(E(graph)$weight)) {
    stop("Edge weights must be non-negative and finite")
  }
  if (length(iter) != 1 || !is.numeric(iter) ||
      !is.finite(iter) || iter != as.integer(iter)) {
    stop("`iter' must be a non-negative finite integer scalar")
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

  ## Project
  graphlets.project(graph, clf, iter)
}

graphlets.project <- function(graph, cliques, iter, Mu=NULL) {

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
  list(Bc=clf[Smb$ix], Muc=Mu[Smb$ix])
}

#################
## Example code

function() {
  library(igraph)
  
  fitandplot <- function(g, gl) {
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

  D1 <- matrix(0, 5, 5)
  D2 <- matrix(0, 5, 5)
  D3 <- matrix(0, 5, 5)
  D1[1:3, 1:3] <- 2
  D2[3:5, 3:5] <- 3
  D3[2:5, 2:5] <- 1
  
  g <- graph.adjacency(D1 + D2 + D3, mode="undirected", weighted=TRUE)
  gl <- graphlets(g, iter=1000)

  fitandplot(g, gl)

  ## Project another graph on the graphlets
  set.seed(42)
  g2 <- set.edge.attribute(g, "weight", value=sample(E(g)$weight))
  gl2 <- graphlets.project(g2, gl$Bc, 1000)
  fitandplot(g2, gl2)
  
}

