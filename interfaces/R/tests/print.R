
library(igraph) ; igraph.options(print.full=TRUE)

g <- graph.ring(5)
summary(g)
g

V(g)$name <- letters[1:vcount(g)]
summary(g)
g

set.seed(42)
E(g)$weight <- runif(ecount(g))
summary(g)
g

g$name <- "A ring"
summary(g)
print(g, v=T)
print(g, e=T)

set.seed(42)
erdos.renyi.game(13, p=0.6, directed=TRUE)

erdos.renyi.game(20, p=0.8)

graph.star(100)

graph.star(100, mode="out")

ba.game(100, m=6, directed=FALSE)

kite <- graph.empty(directed=FALSE) + LETTERS[1:10]
kite <- kite + edges('A','B','A','C','A','D','A','F',
                     'B','D','B','E','B','G', 'C','D','C','F', 
                     'D','E','D','F','D','G', 'E','G', 
                     'F','G','F','H', 'G','H', 'H','I','I','J')
kite
