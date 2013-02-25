
library(igraph)

ce <- simplify(read.graph(gzfile("celegansneural.gml.gz"), format="gml"))
dc <- dyad.census(ce)

round(sum(is.mutual(ce))) == round(dc$mut * 2)
round(ecount(as.undirected(ce, mode="collapse")) - dc$mut) == round(dc$asym)
sum(unlist(dc)) == vcount(ce) * (vcount(ce)-1) / 2

dc
