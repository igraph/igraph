
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

## Eigenvector-centrality, small stress-test

is.principal <- function(M, lambda, eps=1e-12) {
  abs(eigen(M)$values[1] - lambda) < eps
}

is.ev <- function(M, v, lambda, eps=1e-12) {
  max(abs(M %*% v - lambda * v)) < eps
}

is.good <- function(M, v, lambda, eps=1e-12) {
  is.principal(M, lambda, eps) && is.ev(M, v, lambda, eps)
}

for (i in 1:1000) {
  G <- erdos.renyi.game(10, sample(1:20, 1), type="gnm")
  ev <- evcent(G)
  if (!is.good(get.adjacency(G), ev$vector, ev$value)) {
    print("Foobar!")
    break
  }
}

