library(igraph)
g <- \emph{read.graph}("http://localhost/~csardi/karate.net", format="pajek")

community.newman <- function(g) \{
  deg <- \emph{degree}(g) ; ec <- \emph{ecount}(g)
  B <- \emph{get.adjacency}(g) - outer(deg, deg, function(x,y) x*y / 2 / ec)
  diag(B) <- 0
  Re(eigen(B)$vectors[,1])
\}
mem <- community.newman(g)
\emph{V(g)$color} <- ifelse(mem < 0, "grey", "green")

scale <- function(v, a, b) \{ v <- v-min(v) ; v <- v/max(v) ; v <- v * (b-a) ; v+a \}
\emph{V(g)$size} <- scale(abs(mem), 15, 25)
\emph{E(g)$color} <- "grey"
\emph{E(g)[ V(g)[ color=="grey" ] %--% V(g)[ color=="green" ] ]$color} <- "red"
\emph{tkplot}(g, layout=layout.kamada.kawai, vertex.color="a:color",
       vertex.size="a:size", edge.color="a:color")
