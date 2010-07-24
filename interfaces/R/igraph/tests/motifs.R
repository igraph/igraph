
library(igraph)

set.seed(123)
b <- erdos.renyi.game(10000, 4/10000, directed=TRUE)

mno <- graph.motifs.no(b)

mno0 <- graph.motifs.no(b, cut.prob=c(1/3, 0, 0))
mno1 <- graph.motifs.no(b, cut.prob=c(0, 0, 1/3))
mno2 <- graph.motifs.no(b, cut.prob=c(0, 1/3, 0))
c(mno0/mno, mno1/mno, mno2/mno)

mno3 <- graph.motifs.no(b, cut.prob=c(0, 1/3, 1/3))
mno4 <- graph.motifs.no(b, cut.prob=c(1/3, 0, 1/3))
mno5 <- graph.motifs.no(b, cut.prob=c(1/3, 1/3, 0))
c(mno3/mno, mno4/mno, mno5/mno)

######################

set.seed(123)
b <- erdos.renyi.game(10000, 4/10000, directed=TRUE)

m <- graph.motifs(b)

m0 <- graph.motifs(b, cut.prob=c(1/3, 0, 0))
m1 <- graph.motifs(b, cut.prob=c(0, 1/3, 0))
m2 <- graph.motifs(b, cut.prob=c(0, 0, 1/3))
m0/m
m1/m
m2/m

m3 <- graph.motifs(b, cut.prob=c(0, 1/3, 1/3))
m4 <- graph.motifs(b, cut.prob=c(1/3, 1/3, 0))
m5 <- graph.motifs(b, cut.prob=c(1/3, 1/3, 0))
m3/m
m4/m
m5/m

