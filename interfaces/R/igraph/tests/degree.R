
library(igraph)

g <- erdos.renyi.game(100, 1/100)
d <- degree(g)
el <- get.edgelist(g)
all(table(el)==d[d!=0])

g2 <- erdos.renyi.game(100, 2/100, dir=TRUE)
din <- degree(g2, mode="in")
dout <- degree(g2, mode="out")
el2 <- get.edgelist(g2)
all(table(el2[,1]) == dout[dout!=0])
all(table(el2[,2]) == din[din!=0])
