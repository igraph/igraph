
#####################################

R
library(igraph)
library(combinat)

n <- 4
ng <- 2 ^ (n*(n-1))

perms <- permn( 1:n )

to.number <- function(el) {

  g <- graph( t(el), n=n )
  A <- get.adjacency(g)
  diag(A) <- -1
  A <- as.numeric(A)
  A <- A [ A != -1 ]
  sum(A * 2 ^ (0:(length(A)-1)))
  
}

to.edgelist <- function(num, nodes=n) {

  A <- matrix(0, nr=n, nc=n)
  diag(A) <- -1
  v <- numeric()
  
  po <- ng/2
  while (po >= 1) {
    if (num>=po) {
      num <- num - po
      v <- c(v, 1)
    } else {
      v <- c(v, 0)
    }
    po <- po/2
  }
  
  A [ A != -1 ] <- rev(v)
  g <- graph.adjacency(A)
  get.edgelist(g)
}
  

res <- rep(-1, ng)
canon <- numeric()

for (num in 0:(length(res)-1)) {

  print(num)
  el <- to.edgelist(num)
  
  for (p in perms) {
    el2 <- matrix(p[el+1]-1, nc=2)
    num2 <- to.number(el2)
    w <- which(num2 == canon)
    if (length(w) != 0) {
      res[num+1] <- canon[w]
      break
    }
  }
  if (res[num+1]==-1) {
    res[num+1] <- num
    canon <- c(canon, num)
res <- rep(-1, ng)
canon <- numeric()

for (num in 0:(length(res)-1)) {

  print(num)
  el <- to.edgelist(num)
  
  for (p in perms) {
    el2 <- matrix(p[el+1]-1, nc=2)
    num2 <- to.number(el2)
    w <- which(num2 == canon)
    if (length(w) != 0) {
      res[num+1] <- canon[w]
      break
    }
  }
  if (res[num+1]==-1) {
    res[num+1] <- num
    canon <- c(canon, num)
  }
}

  }
}

################################

R
library(igraph)

n <- 3
ng <- 2 ^ (n*(n-1))

to.edgelist <- function(num, nodes=n) {

  A <- matrix(0, nr=n, nc=n)
  diag(A) <- -1
  v <- numeric()
  
  po <- ng/2
  while (po >= 1) {
    if (num>=po) {
      num <- num - po
      v <- c(v, 1)
    } else {
      v <- c(v, 0)
    }
    po <- po/2
  }
  
  A [ A != -1 ] <- rev(v)
  g <- graph.adjacency(A)
  get.edgelist(g)
}

coo <- matrix( c( 0,1, 1,.5, 0,0 ), nc=2, byrow=TRUE)
layout( matrix(1:64, nr=8, nc=8, byrow=TRUE) )
layout.show(64)

for (i in 0:63) {
  
  g <- graph( t(to.edgelist(i)), n=3 )
  par(mar=c(0,0,0,0), mai=c(1,1,1,1)/5)
  class <- graph.isoclass(g)
  col=if (FALSE) "red" else "lightblue"
  plot(g, layout=coo, vertex.size=20, vertex.color=col, labels=NA)
  text(-.9,0,substitute(i %<=>% j, list(i=i, j=class)),
       adj=c(0,.5), cex=2)
}

######################################

coo <- matrix( c(0,1, 1,1, 1,0, 0,0), nc=2, byrow=TRUE)
layout( matrix(1:221, nr=13, nc=17, byrow=TRUE) )
layout.show(221)

for (i in canon) {
  
  g <- graph( t(to.edgelist(i)), n=4 )
  par(mar=c(0,0,0,0), mai=c(1,1,1,1)/5)
  plot(g, layout=coo, vertex.size=40, vertex.color="red", labels=NA,
       edge.color="lightgrey")
  text(0,0,substitute(i, list(i=i)), adj=c(.5,.5), cex=2, col="brown", font=2)
}

##########################################
# UNDIRECTED

R
library(igraph)
library(combinat)

n <- 3
ng <- 2 ^ (n*(n-1)/2)

perms <- permn( 1:n )

to.number <- function(el) {

  g <- graph( t(el), n=n, directed=FALSE )
  A <- get.adjacency(g,type="upper")
  A [ lower.tri(A, diag=TRUE) ] <- -1
  A <- as.numeric(A)
  A <- A [ A != -1 ]
  sum(A * 2 ^ (0:(length(A)-1)))
  
}

to.edgelist <- function(num, nodes=n) {

  A <- matrix(0, nr=n, nc=n)
  A [ lower.tri(A, diag=TRUE) ] <- -1
  v <- numeric()
  
  po <- ng/2
  while (po >= 1) {
    if (num>=po) {
      num <- num - po
      v <- c(v, 1)
    } else {
      v <- c(v, 0)
    }
    po <- po/2
  }
  
  A [ A != -1 ] <- rev(v)
  g <- graph.adjacency(A, mode="upper")
  get.edgelist(g)
}
  

res <- rep(-1, ng)
canon <- numeric()

for (num in 0:(length(res)-1)) {

  print(num)
  el <- to.edgelist(num)
  
  for (p in perms) {
    el2 <- matrix(p[el+1]-1, nc=2)
    num2 <- to.number(el2)
    w <- which(num2 == canon)
    if (length(w) != 0) {
      res[num+1] <- canon[w]
      break
    }
  }
  if (res[num+1]==-1) {
    res[num+1] <- num
    canon <- c(canon, num)
  }
}
