library(igraph)
g <- read.graph("karate.net", format="pajek")

wt <- walktrap.community(g, modularity=TRUE)
dend <- as.dendrogram(wt, use.modularity=TRUE)
plot(dend, nodePar=list(pch=c(NA, 20)))

