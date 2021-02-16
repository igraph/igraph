#######################################################################
## 1) Basics
#######################################################################

## Load the igraph package
library(igraph)

## Create a small graph, A->B, A->C, B->C, C->E, D
## A=0, B=1, C=2, D=3, E=4
g <- graph( c(0,1, 0,2, 1,2, 2,4), n=5 )

## Print a graph to the screen
g

## Create an undirected graph as well
## A--B, A--C, B--C, C--E, D
g2 <- graph( c(0,1, 0,2, 1,2, 2,4), n=5, dir=FALSE )
g2

## How to decide what kind of object a variable refers to
class(g2)
class(1)
class("foobar")

## Is this object an igraph graph?
is.igraph(g)
is.igraph(1:10)

## Summary, number of vertices, edges
summary(g)
vcount(g)
ecount(g)

## Is the graph directed?
is.directed(g)
is.directed(g2)
is.directed(1:10)

## Convert from directed to undirected
as.undirected(g)

## And back
as.directed(as.undirected(g))

## Multiple edges
g <- graph( c(0,1,0,1, 0,2, 1,2, 3,4), n=5 )
g

is.simple(g)
is.multiple(g)

## Remove multiple edges
g <- simplify(g)
is.simple(g)

## Loop edges
g <- graph( c(0,0,0,1, 0,2, 1,2, 3,4), n=5 )
g

is.simple(g)
is.loop(g)

## Remove loop edges
g <- simplify(g)
is.simple(g)

## Naming vertices
g <- graph.ring(10)
V(g)$name <- letters[1:10]
V(g)$name
g
print(g, v=T)

#######################################################################
## 2) Creating graphs
#######################################################################

## The formula interface
library(igraph)

## A simple undirected graph
g <- graph.formula(Alice-Bob-Cecil-Alice, 
                   Daniel-Cecil-Eugene, Cecil-Gordon )

## Another undirected graph, ":" notation
g2 <- graph.formula(Alice-Bob:Cecil:Daniel, 
                    Cecil:Daniel-Eugene:Gordon )

## A directed graph
g3 <- graph.formula(Alice +-+ Bob --+ Cecil 
                    +-- Daniel, Eugene --+ Gordon:Helen )

## A graph with isolate vertices
g4 <- graph.formula(Alice -- Bob -- Daniel, 
                    Cecil:Gordon, Helen )

## "Arrows" can be arbitrarily long
g5 <- graph.formula( Alice +---------+ Bob )

## From edge lists
el <- cbind( c(0, 0, 1, 2),
             c(1, 2, 2, 4) )
g <- graph.edgelist(el)
g

## Symbolic edge lists
el <- cbind( c("Alice", "Alice", "Bob",   "Cecil"),
             c("Bob",   "Cecil", "Cecil", "Ed") )
g <- graph.edgelist(el)
g
summary(g)

## Adjacency matrices
A <- matrix(sample(0:1, 100, rep=TRUE), 10, 10)
g <- graph.adjacency(A)

#######################################################################
## 3) Manipulate graphs, add and delete vertices and edges
#######################################################################

library(igraph)

## Create undirected example graph
g2 <- graph.formula(Alice-Bob:Cecil:Daniel, 
                    Cecil:Daniel-Eugene:Gordon )
print(g2, v=T)

## Remove alice
g3 <- delete.vertices(g2, 0)

## Add three new vertices
g4 <- add.vertices(g3, 3)
print(g4, v=T)
igraph.par("print.vertex.attributes", TRUE)
igraph.par("plot.layout", layout.fruchterman.reingold)
g4

## Add three new vertices, with names this time
g4 <- add.vertices(g3, 3, attr=list(name=c("Helen", "Ike", "Jane")))
g4

## Add some edges as well
g4 <- add.edges(g4, c(6,7, 5,7))
g4

#######################################################################
## 4) Edge sequences
#######################################################################

library(igraph)

## Create directed example graph
g2 <- graph.formula(Alice -+ Bob:Cecil:Daniel, 
                    Cecil:Daniel +-+ Eugene:Gordon )
print(g2, v=T)

## Sequence of all edges
E(g2)

## Edge from a vertex to another
E(g2, P=c(0,1))

## Delete this edge
g3 <- delete.edges(g2, E(g2, P=c(0,1)))
g3

## Get the id of the edge
as.vector(E(g2, P=c(0,1)))

## All adjacent edges of a vertex
E(g2)[ adj(2) ]

## Or multiple vertices
E(g2)[ adj(c(2,0)) ]

## Outgoing edges
E(g2)[ from(2) ]

## Incoming edges
E(g2)[ to(2) ]

## Edges along a path
E(g2, path=c(0,3,4))

#######################################################################
## 5) A real example, creating a graph from data frames
#######################################################################

library(igraph)

## We need some extra code, not in the current version yet
source("http://cneurocvs.rmki.kfki.hu/igraph/plus.R")

## Read the comma-separated value file 
vertices <- read.csv("http://cneurocvs.rmki.kfki.hu/igraph/judicial.csv")

## Read the edges, just a table, space separated 
edges <- read.table("http://cneurocvs.rmki.kfki.hu/igraph/allcites.txt")

## And create the graph
jg <- graph.data.frame(edges, vertices=vertices, dir=TRUE)
summary(jg)

## Save it for future use
save(jg, file="judicial.Rdata.gz", compress=TRUE)

## Easy, huh?

#######################################################################
## 6) Visualizing graphs
#######################################################################

library(igraph)

## Create a scale-free network
g <- barabasi.game(100, m=1)
g <- simplify(g)

## simple plot
igraph.par("plot.layout", layout.fruchterman.reingold)
plot(g, vertex.size=3, vertex.label=NA, edge.arrow.size=0.6)

## interactive
id <- tkplot(g, vertex.size=3, vertex.label=NA, edge.arrow.size=0.6)
coords <- tkplot.getcoords(id)

## 3D
open3d()
rglplot(g, vertex.size=3, vertex.label=NA, edge.arrow.size=0.6)

## A bit better 3D
coords <- layout.kamada.kawai(g, dim=3)
open3d()
rglplot(g, vertex.size=3, vertex.label=NA, edge.arrow.size=0.6, layout=coords)

#######################################################################
## 7) Attributes
#######################################################################

library(igraph)

## Load the jurisdiction network
load("judicial.Rdata.gz")

## If we don't have it then create it again
if (!exists("jg")) {
  source("http://cneurocvs.rmki.kfki.hu/igraph/plus.R")
  vertices <- read.csv("http://cneurocvs.rmki.kfki.hu/igraph/judicial.csv")
  edges <- read.table("http://cneurocvs.rmki.kfki.hu/igraph/allcites.txt")
  jg <- graph.data.frame(edges, vertices=vertices, dir=TRUE)
}

## What do we have?
summary(jg)
V(jg)$year[1:100]
V(jg)$parties[1:10]

## Select vertices based on attributes
V(jg) [ year >= 1990 ]
V(jg) [ overruled!=0 ]

## Group network measures based on attributes
deg.per.year <- tapply(degree(jg, mode="out"), V(jg)$year, mean)

## Plot it
plot( names(deg.per.year), deg.per.year )

## A more advanced example
g <- erdos.renyi.game(100, 1/100)
V(g)$color <- sample( c("red", "black"), 
                     vcount(g), rep=TRUE)
E(g)$color <- "grey"
red <- V(g)[ color == "red" ]
bl <- V(g)[ color == "black" ]
E(g)[ red %--% red ]$color <- "red"
E(g)[ bl  %--% bl ]$color <- "black"
plot(g, vertex.size=5, layout=
     layout.fruchterman.reingold)

#######################################################################
## 8) Centrality and its sensitivity
#######################################################################

library(igraph)

## Create it
g <- graph.formula(Andre----Beverly:Diane:Fernando:Carol,
                   Beverly--Andre:Diane:Garth:Ed,
                   Carol----Andre:Diane:Fernando,
                   Diane----Andre:Carol:Fernando:Garth:Ed:Beverly,
                   Ed-------Beverly:Diane:Garth,
                   Fernando-Carol:Andre:Diane:Garth:Heather,
                   Garth----Ed:Beverly:Diane:Fernando:Heather,
                   Heather--Fernando:Garth:Ike,
                   Ike------Heather:Jane,
                   Jane-----Ike )
g <- simplify(g)

## Hand-drawn coordinates
coords <- c(5,5,119,256,119,256,120,340,478,
            622,116,330,231,116,5,330,451,231,231,231)
coords <- matrix(coords, nc=2)

## Labels the same as names
V(g)$label <- V(g)$name
g$layout <- coords # $

## Take a look at it
plot(g, asp=FALSE, vertex.label.color="blue", vertex.label.cex=1.5,
     vertex.label.font=2, vertex.size=20, vertex.color="white",
     vertex.frame.color="white", edge.color="black")  

## Add degree centrality to labels
V(g)$label <- paste(sep="\n", V(g)$name, degree(g))

## And plot again
plot(g, asp=FALSE, vertex.label.color="blue", vertex.label.cex=1.5,
     vertex.label.font=2, vertex.size=20, vertex.color="white",
     vertex.frame.color="white", edge.color="black")                      

#######################################################################
## 9) Analyzing a (moderately) big graph
#######################################################################

library(igraph)

## Load the jurisdiction network
load("judicial.Rdata.gz")

## If we don't have it then create it again
if (!exists("jg")) {
  source("http://cneurocvs.rmki.kfki.hu/igraph/plus.R")
  vertices <- read.csv("http://cneurocvs.rmki.kfki.hu/igraph/judicial.csv")
  edges <- read.table("http://cneurocvs.rmki.kfki.hu/igraph/allcites.txt")
  jg <- graph.data.frame(edges, vertices=vertices, dir=TRUE)
}

## Basic data
summary(jg)

## Is it a simple graph?
is.simple(jg)

## Is it connected?
is.connected(jg)

## How many components?
no.clusters(jg)

## How big are these?
table(clusters(jg)$csize)

## In-degree distribution
plot(degree.distribution(jg, mode="in"), log="xy")

## Out-degree distribution
plot(degree.distribution(jg, mode="out"), log="xy")

## Largest in- and out-degree, total degree
max(degree(jg, mode="in"))
max(degree(jg, mode="out"))
max(degree(jg, mode="all"))

## Density
graph.density(jg)

## Transitivity
transitivity(jg)

## Transitivity of a random graph of the same size
g <- erdos.renyi.game(vcount(jg), ecount(jg), type="gnm")
transitivity(g)

## Dyad census
dyad.census(jg)

## Triad census
triad.census(jg)

## Authority and Hub scores
authority.score(jg)$vector
cor(authority.score(jg)$vector, V(jg)$auth)

hub.score(jg)$vector
cor(hub.score(jg)$vector, V(jg)$hub)

######################################################################
# 10) Community structure detection
######################################################################

library(igraph)

## The concept of modularity
g <- graph.full(5) %du% graph.full(5)
g$layout <- layout.fruchterman.reingold
plot(g)
V(g)$color <- 2
V(g)[0:4]$color <- 1
plot(g)

## This graph is modular
modularity(g, membership=V(g)$color-1)

## If we have everyone in the same group,
## that is not modular
V(g)$color <- 1
modularity(g, membership=V(g)$color-1)

## If we assign the vertices randomly to two
## groups, that is not very modular, either
V(g)$color <- sample(1:2, vcount(g), rep=TRUE)
plot(g)
modularity(g, membership=V(g)$color-1)

## Keep only the largest component for further analysis
load("judicial.Rdata.gz")

cl <- clusters(jg)
which.max(cl$csize)
cl$membership == which.max(cl$csize)-1
to.keep <- which(cl$membership == which.max(cl$csize)-1)-1
jg2 <- subgraph(jg, to.keep)
summary(jg2)

## Fast & greedy community detection
jg3 <- as.undirected(jg2, mode="collapse")
is.simple(jg3)
## system.time( fc <- fastgreedy.community(jg3) )

## We do it on a smaller graph instead
karate <- read.graph("http://cneurocvs.rmki.kfki.hu/igraph/karate.net",
                     format="pajek")
karate <- simplify(karate)

system.time( fc <- fastgreedy.community(karate) )
fc$modularity
max(fc$modularity)
which.max(fc$modularity)
memb <- community.to.membership(karate, fc$merges, which.max(fc$modularity))
memb

lay <- layout.kamada.kawai(karate)
plot(karate, layout=lay, vertex.size=5, vertex.label=NA,
     vertex.color=memb$membership+1, asp=FALSE)

## Spinglass community detection
## system.time(spc <- spinglass.community(jg3, spins=20))

system.time(spc <- spinglass.community(karate, spins=20))
spc

x11()
plot(karate, layout=lay, vertex.size=5, vertex.label=NA,
     vertex.color=spc$membership+1, asp=FALSE)

######################################################################
# 11) Cohesive blocks
######################################################################

library(igraph)

## Load the graph
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
               12,16, 13,16, 14,15, 15,16,
               17,18,17,19,17,20,
               18,20,18,21,
               19,20,19,22,19,23,
               20,21, 21,22,21,23,
               22,23)-1, dir=FALSE)

## 
V(cb)$label <- seq(vcount(cb))  # $

blocks <- cohesive.blocks(cb)
blocks

summary(blocks)
blocks$blocks
lapply(blocks$blocks, "+", 1)
blocks$block.cohesion
plot(blocks, layout=layout.kamada.kawai,
     vertex.label.cex=2, vertex.size=15,
     vertex.label.color="black")

######################################################################
# 12) Weighted transitivity
######################################################################

library(igraph)

## Define function
wtrans <- function(g) {
  W <- get.adjacency(g, attr="weight")
  WM <- matrix(max(W), nrow(W), ncol(W))
  diag(WM) <- 0
  diag( W %*% W %*% W ) / 
    diag( W %*% WM %*% W)
}

## Try it on an example graph
g <- graph.formula(Andre----Beverly:Diane:Fernando:Carol,
                   Beverly--Andre:Diane:Garth:Ed,
                   Carol----Andre:Diane:Fernando,
                   Diane----Andre:Carol:Fernando:Garth:Ed:Beverly,
                   Ed-------Beverly:Diane:Garth,
                   Fernando-Carol:Andre:Diane:Garth:Heather,
                   Garth----Ed:Beverly:Diane:Fernando:Heather,
                   Heather--Fernando:Garth:Ike,
                   Ike------Heather:Jane,
                   Jane-----Ike )
g <- simplify(g)

## Check the fallback
E(g)$weight <- 1
cbind(transitivity(g, "local"), wtrans(g))

## Now change a weight
E(g)[0]$weight <- 1/2

## Take a look at it
coords <- c(5,5,119,256,119,256,120,340,478,
            622,116,330,231,116,5,330,451,231,231,231)
coords <- matrix(coords, nc=2)
V(g)$label <- paste(sep="", V(g)$name, "\n",
                    round(transitivity(g, "local"),2),
                    ", ", round(wtrans(g), 2))
g$layout <- coords # $
plot(g, asp=FALSE, vertex.label.color="blue", vertex.label.cex=1.5,
     vertex.label.font=2, vertex.size=20, vertex.color="white",
     vertex.frame.color="white", edge.color="black",
     edge.label=E(g)$weight, edge.label.cex=1.5, edge.label.color="red")

######################################################################
# 12) Clique percolation
######################################################################

library(igraph)

## Function to calculate clique communities
clique.community <- function(graph, k) {
  clq <- cliques(graph, min=k, max=k)
  edges <- c()
  for (i in seq(along=clq)) {
    for (j in seq(along=clq)) {
      if ( length(unique(c(clq[[i]], 
             clq[[j]]))) == k+1 ) {
        edges <- c(edges, c(i,j)-1)
      }
    }
  }
  clq.graph <- simplify(graph(edges))
  V(clq.graph)$name <- 
    seq(length=vcount(clq.graph))
  comps <- decompose.graph(clq.graph)
  
  lapply(comps, function(x) {
    unique(unlist(clq[ V(x)$name ]))
  })
}

## Apply it to a graph, this is the example graph from
## the original publication
g <- graph.formula(A-B:F:C:E:D, B-A:D:C:E:F:G, C-A:B:F:E:D, D-A:B:C:F:E,
                   E-D:A:C:B:F:V:W:U, F-H:B:A:C:D:E, G-B:J:K:L:H,
                   H-F:G:I:J:K:L, I-J:L:H, J-I:G:H:L, K-G:H:L:M,
                   L-H:G:I:J:K:M, M-K:L:Q:R:S:P:O:N, N-M:Q:R:P:S:O,
                   O-N:M:P, P-Q:M:N:O:S, Q-M:N:P:V:U:W:R, R-M:N:V:W:Q,
                   S-N:P:M:U:W:T, T-S:V:W:U, U-E:V:Q:S:W:T,
                   V-E:U:W:T:R:Q, W-U:E:V:Q:R:S:T)
g <- simplify(g)

## Hand-made layout to make it look like the original in the paper
lay <- c(387.0763, 306.6947, 354.0305, 421.0153, 483.5344, 512.1145, 
         148.6107, 392.4351, 524.6183, 541.5878, 240.6031, 20, 
         65.54962, 228.0992, 61.9771, 152.1832, 334.3817, 371.8931, 
         421.9084, 265.6107, 106.6336, 57.51145, 605, 20, 124.8780, 
         273.6585, 160.2439, 241.9512, 132.1951, 123.6585, 343.1707, 
         465.1220, 317.561, 216.3415, 226.0976, 343.1707, 306.5854, 
         123.6585, 360.2439, 444.3902, 532.1951, 720, 571.2195, 
         639.5122, 505.3659, 644.3902)
lay <- matrix(lay, nc=2)
lay[,2] <- max(lay[,2])-lay[,2]

## Take a look at it
plot(g, layout=lay, vertex.label=V(g)$name)

## Calculate communities
res <- clique.community(g, k=4)

## Paint them to different colors
colbar <- rainbow( length(res)+1 )
for (i in seq(along=res)) {
  V(g)[ res[[i]] ]$color <- colbar[i+1]
}

## Paint the vertices in multiple communities to red
V(g)[ unlist(res)[ duplicated(unlist(res)) ] ]$color <- "red"

## Plot with the new colors
plot(g, layout=lay, vertex.label=V(g)$name)

######################################################################
# 13) Compare to sna & network
######################################################################

library(igraph)
library(sna)
library(network)

## Some (moderate size) data
load("judicial.Rdata.gz")

## Try using an adjacency matrix
A <- get.adjacency(jg)

## Try a 'network' object
el <- get.edgelist(jg, name=FALSE)
net <- network.initialize(vcount(jg))
net2 <- network::add.edges(net, head=el[,1]+1, tail=el[,2]+1)

## Great! Try calculating degree centrality
deg <- degree(net2)
deg <- igraph::degree(jg)

## Ooops, try eigenvector centrality
ev <- evcent(net2)
ev <- igraph::evcent(jg)$vector

## etc.
