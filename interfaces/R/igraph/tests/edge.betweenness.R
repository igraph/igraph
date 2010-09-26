
library(igraph)

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

bet <- betweenness(kite)
ebet <- edge.betweenness(kite)

bet2 <- sapply(1:vcount(kite), function(x) {
  ae <- E(kite)[ adj(x) ]
  (sum(ebet[ae])-vcount(kite)+1) / 2
})

max(abs(bet - bet2)) < 1e-14

#### Weighted

E(kite)$weight <- sample(1:10, ecount(kite), replace=TRUE)

bet <- betweenness(kite)
ebet <- edge.betweenness(kite)
bet2 <- sapply(1:vcount(kite), function(x) {
  ae <- E(kite)[ adj(x) ]
  (sum(ebet[ae])-vcount(kite)+1) / 2
})

max(abs(bet - bet2)) < 1e-14
