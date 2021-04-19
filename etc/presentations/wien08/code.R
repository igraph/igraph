
## Generating the network figures

## This function creates a lot of graphs

frgraphs <- function() {
  graphs <- list()
  graphs[[1]] <- graph.full(5)
  graphs[[2]] <- graph.full(4)
  graphs[[3]] <- graph.full(6)
  graphs[[4]] <- graph.full(6)
  graphs[[5]] <- graph.full(8)
  graphs[[6]] <- graph.full(2)
  graphs[[7]] <- graph.formula( A-B-C-D-A, A-F, B-G, C-H, D-I )
  graphs[[8]] <- graph.formula( A-B-C-A, D-E-F-D, A-D )
  graphs[[9]] <- graph.ring(5)
  graphs[[10]] <- graph.formula( A-B-C-A, B-D-E-B, C-E-F-C,
                                D-G-H-D, E-H-I-E, F-I-J-F )
  graphs[[11]] <- graph.formula( A-B-C-D-A, A-C, B-D,
                                A-E-F-G-A, A-F, E-G,
                                B-H-I-J-B, B-I, H-J,
                                C-K-L-M-C, C-L, K-M,
                                D-N-O-P-D, D-O, N-P )
  graphs[[12]] <- graph.formula( A-B-C-D-F-E-A, A-D, E-C, B-F)
  graphs[[13]] <- graph.formula( A-B:D:H, B-F:G, C-D:E:I, D-F, E-F:L, F-G:I:J:M,
                                G-N, H-K:L, I-K, J-M:N, K-M, L-N )
  graphs[[14]] <- simplify(graph.formula( A-B:H:L, B-A:C:K, C-B:D:M, D-C:E:L:M:N,
                                         E-D:F:G, F-E:I:L, G-E:H:N, H-A:G:J,
                                         I-F:J:K, J-H:I:M, K-B:I:N, L-A:D:F,
                                         M-C:D:J, N-D:G:K ))
  graphs[[15]] <- graph.formula( A-B-C-D-E-F-G-H-I-A,
                                J-E:A:F, K-A:F:B, L-F:B:G, M-B:G:C, N-G:C:H,
                                O-C:H:D, P-H:D:I, Q-D:I:E, R-I:E:A )
  graphs[[16]] <- graph.formula( 0-1:4:8:11:15:19:22, 1-0:2:6:24, 2-1:3:5:23:29:30,
                                3-2, 4-0:5, 5-2:4:6:9, 6-1:5:8:31:7:30, 7-6:30:31,
                                8-0:6:9:13, 9-5:8:12:32:10:31, 10-9, 11-0:12,
                                12-9:11:26:13, 13-8:12:27:33:14:32, 14-13:32:33,
                                15-0:16, 16-15:20:17:26, 17-16:19:35:18:34:27,
                                18-17:34:35, 19-0:24:20:17, 20-19:23:36:21:35:16,
                                21-20, 22-0:23, 23-22:2:24:20, 24-23:1:29:25:36:19,
                                25-24:29:36, 26-16:12:33:27:34, 27-26:13:28:17,
                                28-27, 29-2:24:25, 30-2:6:7, 31-6:7:9,
                                32-9:13:14, 33-13:14:26, 34-26:17:18, 35-17:18:20,
                                36-20:24:25 )
  graphs[[16]] <- simplify(graphs[[16]])
  return graphs
}

#####################

library(igraph)
library(Cairo)

graphs <- frgraphs()
lay <- lapply(graphs, layout.fruchterman.reingold, niter=3000)

CairoPDF(file="frplots.pdf")
par(mai=c(0,0,0,0))
layout(matrix(1:16, nr=4, byrow=TRUE))
for (i in seq(along=graphs)) {
  plot(graphs[[i]], layout=lay[[i]],
       vertex.label=NA, vertex.size=13, edge.color="black",
       vertex.color="red")
}
dev.off()

#####################

library(igraph)
library(Cairo)

g <- erdos.renyi.game(300, 250, type="gnm")
colors <- rainbow( length(clusters(g)$csize) )
V(g)$color <- colors[ clusters(g)$membership+1 ]
CairoPDF(file="ercomps.pdf")
par(mai=c(0,0,0,0))
plot(g, layout=layout.fruchterman.reingold, vertex.label=NA,
     vertex.size=4, edge.color="black")
dev.off()

######################

library(igraph)
library(rgl)

g <- graph.lattice( c(5,5,5), directed=TRUE, mutual=TRUE )
V(g)$color <- "red"
l <- layout.kamada.kawai(g, dim=3)
open3d()
rglplot(g, vertex.label=NA, layout=l)
Sys.sleep(2) ; rgl.snapshot(file="3dplot.png")

######################
# Draw the example graph

g <- graph.formula( A-B:C, B-C, C-E, D )
V(g)$label <- V(g)$name
CairoPDF(file="small1.pdf")
plot(g, layout=layout.fruchterman.reingold, frame=TRUE,
     vertex.size=50, vertex.label.cex=3, edge.color="black")
dev.off()

CairoPDF(file="small2.pdf")
plot(g, layout=layout.fruchterman.reingold, frame=TRUE,
     vertex.size=50, vertex.label.cex=3, edge.color="black")
dev.off()

l <- layout.fruchterman.reingold(g)

CairoPDF(file="small3.pdf")
plot(g, layout=l, frame=TRUE,
     vertex.size=50, vertex.label.cex=3, edge.color="black")
dev.off()

g2 <- graph.formula( A-+B:C, B-+C, C-+E, D )
V(g2)$label <- V(g)$name
CairoPDF(file="small4.pdf")
plot(g2, layout=l, frame=TRUE,
     vertex.size=50, vertex.label.cex=3, edge.color="black")
dev.off()

#####################

# graph representation

CairoPDF(file="adjacency.pdf")
par(mai=c(0,0,0,0))
plot(NA, type="n", xlim=c(-2,10), ylim=c(0,12), frame=FALSE)

for (i in seq(0.5,9.5,by=1)) {
  symbols(x=rep(i,10), y=seq(0.5,9.5,by=1), squares=rep(1,10),
          inches=FALSE, add=TRUE)
}

names <- c("Alice", "Bob", "Cecil", "Diana", "Esmeralda",
           "Fabien", "Gigi", "Helen", "Iannis", "Jennifer")

text(0, seq(9.5,0.5,by=-1), pos=2, names)
text( seq(0.3,9.3,by=1), 10.2, pos=4, names, srt=90)

g <- erdos.renyi.game(10, 20, type="gnm")
A <- get.adjacency(g)

x <- seq(0.5,9.5,by=1)
for (i in 1:10) {
  text(x=rep(x[i],10), y=rev(x), A[i,], cex=1.5)
}
dev.off()

##########

V(g)$name <- names
el <- get.edgelist(g)
# el <- el[1:16,]
# el[16,] <- c("etc.", "")

CairoPDF(file="edgelist.pdf", height=5, width=2)
par(mai=c(0,0,0,0))
plot(NA, type="n", xlim=c(0,4), ylim=c(0,nrow(el)), frame=FALSE)

text(rep(0,nrow(el)), seq(nrow(el)-0.5, 0.5, by=-1), el[,1], pos=4)
text(rep(2,nrow(el)), seq(nrow(el)-0.5, 0.5, by=-1), el[,2], pos=4)
dev.off()

CairoPNG(file="edgelist.png", height=500, width=200)
par(mai=c(0,0,0,0))
plot(NA, type="n", xlim=c(0,4), ylim=c(0,nrow(el)), frame=FALSE)

text(rep(0,nrow(el)), seq(nrow(el)-0.5, 0.5, by=-1), el[,1], pos=4)
text(rep(2,nrow(el)), seq(nrow(el)-0.5, 0.5, by=-1), el[,2], pos=4)
dev.off()

############

adj <- get.adjlist(g)

# x11(height=5, width=6)
CairoPDF(file="adjlist.pdf", height=5, width=6)
par(mai=c(0,0,0,0))
plot(NA, type="n", xlim=c(0,10), ylim=c(0,length(adj)), frame=FALSE)

text(rep(0,length(adj)), seq(length(adj)-0.5, 0.5, by=-1), V(g)$name, pos=4)
for (i in seq(adj)) {
  tt <- paste(V(g)$name[ adj[[i]]+1 ], collapse=", ")
  text(rep(2, length(adj[[i]])), rev(x)[i], tt, pos=4)
}
dev.off()

################

l <- layout.kamada.kawai(g)

CairoPDF(file="example.pdf")
par(mai=c(0,0,0,0))
plot(g, layout=l, vertex.size=30,
     vertex.label=V(g)$name, vertex.label.cex=1.5,
     vertex.color="white", vertex.frame.color="white")
dev.off()

CairoPDF(file="ex-deg.pdf")
par(mai=c(0,0,0,0))
plot(g, layout=l, vertex.size=30,
     vertex.label=paste(V(g)$name, degree(g), sep="\n"),
     vertex.label.cex=1.5,
     vertex.color="white", vertex.frame.color="white")
dev.off()

CairoPDF(file="ex-close.pdf")
par(mai=c(0,0,0,0))
plot(g, layout=l, vertex.size=30,
     vertex.label=paste(V(g)$name, round(closeness(g),2), sep="\n"),
     vertex.label.cex=1.5,
     vertex.color="white", vertex.frame.color="white")
dev.off()

CairoPDF(file="ex-betw.pdf")
par(mai=c(0,0,0,0))
plot(g, layout=l, vertex.size=30,
     vertex.label=paste(V(g)$name, round(betweenness(g),2), sep="\n"),
     vertex.label.cex=1.5,
     vertex.color="white", vertex.frame.color="white")
dev.off()

CairoPDF(file="ex-ev.pdf")
par(mai=c(0,0,0,0))
plot(g, layout=l, vertex.size=30,
     vertex.label=paste(V(g)$name, round(evcent(g)$vector, 2), sep="\n"),
     vertex.label.cex=1.5,
     vertex.color="white", vertex.frame.color="white")
dev.off()

pg <- page.rank(g)$vector
pg <- round(pg/max(pg), 2)
CairoPDF(file="ex-pagerank.pdf")
par(mai=c(0,0,0,0))
plot(g, layout=l, vertex.size=30,
     vertex.label=paste(V(g)$name, pg, sep="\n"),
     vertex.label.cex=1.5,
     vertex.color="white", vertex.frame.color="white")
dev.off()

#########################
## Cohesive blocking

cb <- graph( c(1,2,1,3,1,4,1,5,1,6,
              2,3,2,4,2,5,2,7,
              3,4,3,6,3,7,
              4,5,4,6,4,7,
              5,6,5,7,5,21,
              6,7,
              7,8,7,11,7,14,7,19,
              8,9,8,11,8,14,
              9,10,
              10,12,10,13,
              11,12,11,14,
              12,16,
              13,16,
              14,15,
              15,16,
              17,18,17,19,17,20,
              18,20,18,21,
              19,20,19,22,19,23,
              20,21,
              21,22,21,23,
              22,23)-1, dir=FALSE)

V(cb)$label <- seq(vcount(cb))

