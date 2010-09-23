
library(igraph)

## Dense

I <- matrix(sample(0:1, 35, replace=TRUE, prob=c(3,1)), nc=5)
g <- graph.incidence(I)
I2 <- get.incidence(g)
all(I==I2)
all(rownames(I2) == 1:7)
all(colnames(I2) == 8:12)

## Sparse

I3 <- get.incidence(g, sparse=TRUE)
all(I==I3)
all(rownames(I3) == 1:7)
all(colnames(I3) == 8:12)
