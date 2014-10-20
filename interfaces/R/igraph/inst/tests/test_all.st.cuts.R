
context("all.st.cuts")

test_that("all.st.cuts works", {

  library(igraph)

  unvs <- function(x) lapply(x, as.vector)

  g <- graph_from_literal( a -+ b -+ c -+ d -+ e )
  cc <- st_cuts(g, source="a", target="e")
  expect_that(unvs(cc$cuts), equals(list(1,2,3,4)))
  expect_that(unvs(cc$partition1s), equals(list(1, 1:2, 1:3, 1:4)))

  g2 <- graph_from_literal( s -+ a:b -+ t, a -+ 1:2:3 -+ b )
  cc <- st_cuts(g2, source="s", target="t")
  expect_that(unvs(cc$cuts), equals(list(c(1,2), c(1,7), c(2,3,4,5,6),
                                   c(2,3,4,5,10), c(2,3,4,6,9),
                                   c(2,3,4,9,10), c(2,3,5,6,8),
                                   c(2,3,5,8,10), c(2,3,6,8,9),
                                   c(2,3,8,9,10), c(3,7))))
  expect_that(unvs(cc$partition1s),
              equals(list(1, c(1,3), c(1,2), c(1,2,7), c(1,2,6),
                          c(1,2,6,7), c(1,2,5), c(1,2,5,7), c(1,2,5,6),
                          c(1,2,5,6,7), c(1,2,5,6,7,3))))

  g3 <- graph_from_literal( s -+ a:b -+ t, a -+ 1:2:3:4:5 -+ b )
  cc <- st_min_cuts(g2, source="s", target="t")
  expect_that(cc$value, equals(2))
  expect_that(unvs(cc$cuts), equals(list(c(1,2), c(1,7), c(3,7))))
  expect_that(unvs(cc$partition1s), equals(list(1, c(1,3), c(1,3,2,7,6,5))))
})
