
R

library(igraph)

# KARATE network, we will need this

karate <- graph( c( 0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
                   0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
                   0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
                   0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
                   1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
                   2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
                   2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
                   4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
                   6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
                   13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
                   18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
                   22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
                   23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
                   25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
                   28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
                   31, 32, 31, 33, 32, 33), directed=FALSE)

# DIAMETER OF THE KARATE CLUB

g <- karate
d <- get.diameter(g)
E(g)$color <- "grey"
E(g)$width <- 1
E(g, path=d)$color <- "red"
E(g, path=d)$width <- 2
V(g)$label.color <- "blue"
V(g)$color  <- "SkyBlue2"
V(g)[ d ]$label.color <- "black"
V(g)[ d ]$color <- "red"
plot(g, layout=layout.fruchterman.reingold, 
     vertex.label.dist=0, vertex.size=15)
title(main="Diameter of the Zachary Karate Club network",
      xlab="created by igraph 0.4")
axis(1, labels=FALSE, tick=TRUE)
axis(2, labels=FALSE, tick=TRUE)


# DIAMETER OF A SCALE-FREE GRAPH

g <- barabasi.game(100, directed=FALSE)
d <- get.diameter(g)
E(g)$color <- "SkyBlue2"
E(g)$width <- 1
E(g, path=d)$color <- "red"
E(g, path=d)$width <- 2
V(g)$label.color <- V(g)$color  <- "blue"
V(g)[ d ]$label.color <- V(g)[ d ]$color <- "red"
plot(g, layout=layout.fruchterman.reingold, 
     vertex.label.dist=0.6, vertex.size=3)
title(main="Diameter of a small scale-free random graph",
      xlab="created by igraph 0.4")
axis(1, labels=FALSE, tick=TRUE)
axis(2, labels=FALSE, tick=TRUE)

# TKPLOT, NEWMAN-CLUSTERED KARATE-CLUB

g <- karate

community.newman <- function(g) {
  d <- degree(g)
  B <- get.adjacency(g)-outer(d, d, function(x,y) x*y/2/ecount(g))
  diag(B) <- 0
  eigen(B)$vectors[,1]
}

mem <- community.newman(g)
V(g)$color <- ifelse(mem < 0, "blue", "green")
V(g)$size <- abs(mem) * 35
E(g)$color <- "darkblue"
E(g) [ V(g)[ color=="green" ] %--% V(g)[ color=="green" ] ]$color <- "darkolivegreen"
E(g) [ V(g)[ color=="blue" ] %--% V(g)[ color=="green" ] ]$color <- "red"
tkplot(g, layout=layout.fruchterman.reingold,
       vertex.label.dist=1)

# A 3D PLOT

g <- karate
l <- layout.kamada.kawai(g, dim=3)
wt <- walktrap.community(g)
mem <- community.to.membership(g, wt$merges, steps=30)$membership
cb <- heat.colors(4)
E(g)$color <- "black"
for (i in 0:3) {
  V(g)[mem==i]$color <- cb[i+1]
  E(g)[ V(g)[mem==i] %--% V(g)[mem==i] ]$color <- cb[i+1]
}

open3d()
rglplot(g, layout=l, vertex.label.dist=1)
rgl.postscript("karate3d.eps")

# DENDROGRAM, KARATE NETWORK

g <- karate

wt <- walktrap.community(g, modularity=TRUE)
dend <- as.dendrogram(wt, use.modularity=TRUE)
xfig(file="dendrogram.fig")
plot(dend, nodePar=list(pch=c(NA, 20)))
dev.off()

# SOME DEGREE DISTRIBUTIONS

g <- barabasi.game(100000)
d <- degree(g, mode="in")
dd <- degree.distribution(g, mode="in", cumulative=TRUE)
alpha <- power.law.fit(d, xmin=20)
plot(dd, log="xy", xlab="degree", ylab="cumulative frequency",
     col="red", main="Nonlinear preferential attachment")
lines(10:500, 10*(10:500)^(-coef(alpha)+1))

g <- barabasi.game(100000, power=0.9)
dd <- degree.distribution(g, mode="in", cumulative=TRUE)
points(dd, col="green", pch=2)

g <- barabasi.game(100000, power=0.8)
dd <- degree.distribution(g, mode="in", cumulative=TRUE)
points(dd, col="blue", pch=3)

g <- barabasi.game(100000, power=0.7)
dd <- degree.distribution(g, mode="in", cumulative=TRUE)
points(dd, col="orange", pch=4)

g <- barabasi.game(100000, power=0.6)
dd <- degree.distribution(g, mode="in", cumulative=TRUE)
points(dd, col="purple", pch=5)

legend(1, 1e-5, c("1.0", "0.9", "0.8", "0.7", "0.6"),
       col=c("red", "green", "blue", "orange", "purple"),
       pch=1:5, ncol=1, yjust=0, lty=0)

       
