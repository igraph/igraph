
context("closeness")

test_that("closeness works", {
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

  clo <- closeness(kite) * (vcount(kite)-1)
  expect_that(round(sort(clo, decreasing=TRUE), 3),
              equals(c(Fernando=0.643, Garth=0.643, Diane=0.600,
                       Heather=0.600, Andre=0.529, Beverly=0.529,
                       Carol=0.500, Ed=0.500, Ike=0.429, Jane=0.310)))

  clo2 <- closeness(kite, normalized=TRUE)
  expect_that(clo, equals(clo2))
})

## TODO: weighted closeness

test_that("closeness centralization works", {

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

  c1 <- closeness(kite, normalized=TRUE)
  c2 <- centralization.closeness(kite)
  expect_that(unname(c1), equals(c2$res))
  expect_that(c2$centralization, equals(0.270374931581828))
  expect_that(c2$theoretical_max, equals(4.23529411764706))
  
})
