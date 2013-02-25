
library(igraph)

gc <- function(graph) {
  clu <- clusters(graph)
  induced.subgraph(graph, which(clu$membership==which.max(clu$csize)))
}

g <- gc(erdos.renyi.game(30, 8/30))
ec <- edge.connectivity(g)
ecST <- Inf
for (j in 1:(vcount(g)-1)) {
  for (k in (j+1):vcount(g)) {
    ec2 <- edge.connectivity(g, source=j, target=k)
    if (ec2 < ecST) { ecST <- ec2 }
  } 
}
ec == ecST

####

kite <- graph.formula(Andre    - Beverly:Carol:Diane:Fernando,
                      Beverly  - Andre:Diane:Ed:Garth,
                      Carol    - Andre:Diane:Fernando,
                      Diane    - Andre:Beverly:Carol:Ed:Fernando:Garth,
                      Ed       - Beverly:Diane:Garth,
                      Fernando - Andre:Carol:Diane:Garth:Heather,
                      Garth    - Beverly:Diane:Ed:Fernando:Heather,
                      Heather  - Fernando:Garth:Ike,
                      Ike      - Heather:Jane,
                      Jane     - Ike)
kite <- simplify(kite)
edge.connectivity(kite, source="Heather", target="Andre")
edge.connectivity(kite, source="Garth", target="Andre")
edge.connectivity(kite, source="Garth", target="Ike")

