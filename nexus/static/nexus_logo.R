
library(igraph)

g <- graph.star(5, mode="undirected")
g$layout <- rbind( c(0,0), c(-1,1), c(1,1), c(1,-1), c(-1,-1) )
V(g)$size <- 80
E(g)$width <- 30
E(g)$color <- "black"
V(g)$label <- c("X", "N", "E", "S", "U")
V(g)$color <- "gold"
V(g)$label.cex <- 7

pdf("nexus_logo.pdf")
par(mar=c(0,0,0,0), bg="lightgrey")
plot(g)
dev.off()

pdf("nexus_logo_small.pdf")
par(mar=c(0,0,0,0), bg="lightgrey")
plot(g, vertex.label=NA)
dev.off()

pdf("nexus_logo_black.pdf")
par(mar=c(0,0,0,0), bg="lightgrey")
plot(g, vertex.color="black", vertex.label.color="white")
dev.off()

pdf("nexus_logo_black_small.pdf")
par(mar=c(0,0,0,0), bg=rgb(0.9,0.9,0.9))
plot(g, vertex.color="black", vertex.label.color="white", vertex.label=NA)
dev.off()

png("nexus_logo_black_small.png", width=160, height=160)
par(mar=c(0,0,0,0), bg="transparent")
plot(g, vertex.color="black", vertex.label.color="white", vertex.label=NA)
dev.off()
