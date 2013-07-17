
context("graph.adhesion")

test_that("graph.adhesion works", {

  library(igraph)

  g <- graph.famous("Zachary")
  expect_that(graph.adhesion(g), equals(1))
  expect_that(graph.cohesion(g), equals(1))

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

  expect_that(graph.adhesion(kite), equals(1)) 
  expect_that(graph.cohesion(kite), equals(1))

  camp <- graph.formula(Harry:Steve:Don:Bert - Harry:Steve:Don:Bert,
                        Pam:Brazey:Carol:Pat - Pam:Brazey:Carol:Pat,
                        Holly   - Carol:Pat:Pam:Jennie:Bill,
                        Bill    - Pauline:Michael:Lee:Holly,
                        Pauline - Bill:Jennie:Ann,
                        Jennie  - Holly:Michael:Lee:Ann:Pauline,
                        Michael - Bill:Jennie:Ann:Lee:John,
                        Ann     - Michael:Jennie:Pauline,
                        Lee     - Michael:Bill:Jennie,
                        Gery    - Pat:Steve:Russ:John,
                        Russ    - Steve:Bert:Gery:John,
                        John    - Gery:Russ:Michael)

  expect_that(graph.adhesion(camp), equals(2))
  expect_that(graph.cohesion(camp), equals(2))
})
