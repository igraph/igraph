
library(igraph)

g <- erdos.renyi.game(100, 1/100)
d <- degree(g)
el <- get.edgelist(g)
all(table(el)==d[d!=0])

all(degree(g) / (vcount(g)-1) == degree(g, normalized=TRUE))

g2 <- erdos.renyi.game(100, 2/100, dir=TRUE)
din <- degree(g2, mode="in")
dout <- degree(g2, mode="out")
el2 <- get.edgelist(g2)
all(table(el2[,1]) == dout[dout!=0])
all(table(el2[,2]) == din[din!=0])

all(degree(g2, mode="in") / (vcount(g2)-1) ==
    degree(g2, mode="in", normalized=TRUE))
all(degree(g2, mode="out") / (vcount(g2)-1) ==
    degree(g2, mode="out", normalized=TRUE))
all(degree(g2, mode="all") / (vcount(g2)-1) ==
    degree(g2, mode="all", normalized=TRUE))
