library(igraph)

g <- read.graph("karate.net", format="pajek")
d <- get.diameter(g)

E(g)$color <- "grey"
E(g)$width <- 1
E(g, path=d)$color <- "red"
E(g, path=d)$width <- 2
V(g)$label.color <- "blue"
V(g)$color  <- "SkyBlue2"
V(g)[ d ]$label.color <- "black"
V(g)[ d ]$color <- "red"

plot(g, layout=layout.fruchterman.reingold, 
     vertex.label.dist=0, vertex.size=15)
title(main="Diameter of the Zachary Karate Club network",
      xlab="created by igraph 0.4")
axis(1, labels=FALSE, tick=TRUE)
axis(2, labels=FALSE, tick=TRUE)

