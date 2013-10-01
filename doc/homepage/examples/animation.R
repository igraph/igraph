library(igraph)

g <- read.graph("karate.net", format="pajek")
l <- layout.kamada.kawai(g, niter=1000)
ebc <- edge.betweenness.community(g)

colbar <- rainbow(6)
colbar2 <- c(rainbow(5), rep("black",15))

for (i in 1:20) {

  g2 <- delete.edges(g, ebc$removed.edges[seq(length=i-1)])
  
  eb <- edge.betweenness(g2)
  cl <- clusters(g2)$membership
  q <- modularity(g, cl)

  E(g2)$color <- "grey"
  E(g2)[ order(eb, decreasing=TRUE)[1:5] ]$color <- colbar2[1:5]

  E(g2)$width <- 1
  E(g2)[ color != "grey" ]$width <- 2

  # png(sprintf("eb-community-%04d.png", i))
  plot(g2, layout=l, vertex.size=6, vertex.label=NA,
       edge.label.color="red", vertex.color=colbar[cl+2],
       edge.label.font=2)
  title(main=paste("Q=", round(q,3)), font=2)
  ty <- seq(1,by=-strheight("1")*1.5, length=20)
  text(-1.3, ty, adj=c(0,0.5), round(sort(eb, dec=TRUE)[1:20],2), col=colbar2,
       font=2)

  # dev.off()
  if (interactive()) Sys.sleep(1)
}
