
library(igraph)

set.seed(42)

g <- barabasi.game(10, m=3, algorithm="bag")
is.multiple(g)
count.multiple(g)
is.multiple(simplify(g))
all(count.multiple(simplify(g)) == 1)
     
## Direction of the edge is important
is.multiple(graph( c(1,2, 2,1) ))
is.multiple(graph( c(1,2, 2,1), dir=FALSE ))

## Remove multiple edges but keep multiplicity
g <- barabasi.game(10, m=3, algorithm="bag")
E(g)$weight <- 1
g <- simplify(g)
any(is.multiple(g))
E(g)$weight

