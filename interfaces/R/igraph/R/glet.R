
## Algorithm 1:
## D is NxN non-negative weighted matrix
## level is a number between minimum and maximum element in D
## puts a Threshold on D at level and returns all maximal cliques of
## the binary network 

threshold.net <- function(D, level) {
  N <- nrow(D)
  D.t <- data.matrix(D >= level)
  El <- which(D.t > 0, arr.ind=TRUE)
  Dth <- graph.edgelist(El, directed=FALSE)

  clqt <- maximal.cliques(Dth)
  pobt <- length(clqt)

  Bc <- matrix(0, N, pobt)
  for (j in 1:pobt) { Bc[clqt[[j]], j] <- 1 }

  W <- matrix(0, 1, pobt)
  for (j in 1:pobt) { W[j] <- sum(Bc[, j]) }
  WS <- sort(W, decreasing = TRUE, index.return = TRUE)
  Bc <- Bc[, WS$ix[1:length(WS$ix)]]
  if (pobt == 1) { Bc <- t(Bc) }

  Bc
}


## Algorithm 2:
## D is NxN non negative weighted matrix
##  projects the network on the basis elements(Bc) and finds Mu

project.net <- function(D, Bc, BB, Mu, iter) {
  K <- ncol(Bc)
  Qt <- list
  a <- 0
  for (j in 1:K) { a[j] <- sum(Bc[, j])^2 }

  for (i in 1:iter) {
    Dh <- Bc %*% diag(Mu) %*% t(Bc)
    for (j in 1:K) {
      Qt <- BB[which(BB[, 2] == j), 1]
      Mu[j] <- Mu[j] * sum(D[Qt, Qt] / (Dh[Qt, Qt] + .0001)) / a[j]
    }
  }
  Mu
}


## Main:
## Finds all basis(Bc) for the network D (using algorithm 1), projects
## on the Bc and finds Muc (using algorithm 2). Input D must be
## symettric with all nonnegative elements

graphlets <- function(D, iter) {  

  print("Summary of the network:")
  Summary.net(D)
  print("Bit String Decomposition:")
  N <- nrow(D)
  Bc <- matrix(0, N, 1)
  ls <- sort(D, decreasing=TRUE)
  ldif <- diff(ls) < 0
  L <- ls[which(ldif)]
  L <- c(L, 0.0001)

  for (level in L) {
    Bt <- threshold.net(D, level)
    pobt <- ncol(Bt)
    W <- matrix(0, 1, pobt)
    for (j in 1:pobt) { W[j] <- sum(Bt[, j]) }
    WS <- sort(W, decreasing=TRUE, index.return=TRUE)
    Bt <- Bt[, WS$ix[1:length(WS$ix)]]
    ## Bt <- Bt[, which(WS$x > 1)]

    Wc <- matrix(0, 1, ncol(Bc))
    for (j in 1:ncol(Bc)) { Wc[j] <- sum(Bc[, j] * (1:N)) }

    indic <- length(Bt) / N

    if (indic > 1) {
      for (j in 1:ncol(Bt)) {
        flag <- 0
        for (iij in which(Wc == sum(Bt[, j] * (1:N)))) {
          if (sum(Bc[, iij] == Bt[, j]) == N) { flag <- 1 }
        }
        if ((flag == 0) & (sum(Bt[, j]) > 1) ) {
          Bc <- matrix(c(Bc, Bt[, j]), nrow=N)
        }
      }
    }
    
  }
  
  Bc <- Bc[, 2:ncol(Bc)]
  BB <- which(Bc == 1, arr.ind=TRUE)
  Mu <- array(1, ncol(Bc))
  ## iter=10000
  Mu <- project.net(D, Bc, BB, Mu, iter)
  Smb <- sort(Mu, decreasing=TRUE, index=TRUE)

  list(Bc=Bc[, Smb$ix], Muc=Mu[Smb$ix])
}

Summary.net <- function(Data.adj) {
  ## print("Summary of network")
  out=sprintf("Node= %d", nrow(Data.adj))
  print(out)
  out=sprintf("Edge= %f", sum(Data.adj > 0) / 2)
  print(out)
  ## out=sprintf("Toal msg= %f",sum(Data.adj)/2)
  ## print(out)
  out=sprintf("average weight per edge= %f", sum(Data.adj) / sum(Data.adj > 0))
  print(out)
  out=sprintf("max weight= %f", max(Data.adj))
  print(out)
  ## pdf("/n/airoldi_lab/azari/MsgFb/Results/Alg1/Weights.pdf")      
  ## hist(Data.adj[which(Data.adj>0)])
  ## dev.off() 
}

#################
## Example code

function() {
  D1 <- matrix(0, 5, 5)
  D2 <- matrix(0, 5, 5)
  D3 <- matrix(0, 5, 5)
  D1[1:3, 1:3] <- 2
  D2[3:5, 3:5] <- 3
  D3[2:5, 2:5] <- 1
  
  D <- D1 + D2 + D3
  
  gl <- graphlets(D, iter=1000)
  
  g <- graph.adjacency(D, weighted=TRUE, mode="undirected", diag=FALSE)
  V(g)$color <- "white"
  E(g)$label <- E(g)$weight
  E(g)$label.cex <- 2
  E(g)$color <- "black"

  plot.new()
  layout(matrix(1:6, nrow=2, byrow=TRUE))
  co <- layout.kamada.kawai(g)
  par(mar=c(1,1,1,1))
  plot(g, layout=co)
  for (i in 1:ncol(gl$Bc)) {
    sel <- gl$Bc[,i] != 0
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

