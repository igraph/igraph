
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
round(evcent(kite)$vector, 3)

