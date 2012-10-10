
library(igraph) 

########################################
## Disjoint union, named graphs


g1 <- g2 <- graph.ring(10)
g1$foo <- "bar"
V(g1)$name <- letters[ 1:10]
V(g2)$name <- letters[11:20]
E(g1)$weight <- 1:10
E(g2)$weight <- 10:1

V(g1)$a1 <- 1:10
V(g2)$a2 <- 11:20

E(g1)$b1 <- 1:10
E(g2)$b2 <- 11:20

g1 + g2

V(g1+g2)$name
V(g1+g2)$a1
V(g1+g2)$a2

E(g1+g2)$weight
E(g1+g2)$b1
E(g1+g2)$b2

########################################
## Union, unnamed graphs

library(igraph)

g1 <- graph.ring(10)
g2 <- graph.ring(13)
g1$foo <- "bar"
E(g1)$weight <- 1:10
E(g2)$weight <- 13:1

V(g1)$a1 <- 1:10
V(g2)$a2 <- 11:23

E(g1)$b1 <- letters[1:10]
E(g2)$b2 <- letters[11:23]

g <- graph.union(g1, g2)
g

V(g)$a1
V(g)$a2

E(g)$b1
E(g)$b2

df1 <- get.data.frame(g)
df1 <- df1[ order(df1$from, df1$to), c(1,2,3,5,4,6)]
df2 <- merge(get.data.frame(g1), get.data.frame(g2),
             by=c("from", "to"), all=TRUE)
all(is.na(df1) == is.na(df2)) && all(df1 == df2, na.rm=TRUE)

###########################################
## Union named graphs

library(igraph)

g1 <- graph.ring(10)
g2 <- graph.ring(13)
V(g1)$name <- letters[seq_len(vcount(g1))]
V(g2)$name <- letters[seq_len(vcount(g2))]

g1$foo <- "bar"
E(g1)$weight <- 1:10
E(g2)$weight <- 13:1

V(g1)$a1 <- 1:10
V(g2)$a2 <- 11:23

E(g1)$b1 <- letters[1:10]
E(g2)$b2 <- letters[11:23]

g <- graph.union(g1, g2)
g
V(g)$name

gg <- graph.union(g1, g2, keep.all.vertices=FALSE)
gg
V(gg)$name

get.data.frame(g, what="both")
get.data.frame(gg, what="both")
