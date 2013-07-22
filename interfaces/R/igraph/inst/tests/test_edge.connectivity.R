
context("edge.connectivity")

test_that("edge.connectivity works", {

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
  expect_that(ec, equals(ecST))

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

  ec1 <- edge.connectivity(kite, source="Heather", target="Andre")
  ec2 <- edge.connectivity(kite, source="Garth", target="Andre")
  ec3 <- edge.connectivity(kite, source="Garth", target="Ike")
  expect_that(ec1, equals(2))
  expect_that(ec2, equals(4))
  expect_that(ec3, equals(1))
})
