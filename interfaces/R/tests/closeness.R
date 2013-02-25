
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
clo <- structure(closeness(kite) * (vcount(kite)-1), names=V(kite)$name)
round(sort(clo, decreasing=TRUE), 3)

clo2 <- closeness(kite, normalized=TRUE)
all(clo==clo2)

## TODO: weighted closeness
