
context("all.st.cuts")

test_that("all.st.cuts works", {

  library(igraph)

  g <- graph.formula( a -+ b -+ c -+ d -+ e )
  cc <- stCuts(g, source="a", target="e")
  expect_that(cc$cuts, equals(list(1,2,3,4)))
  expect_that(cc$partition1s, equals(list(1, 1:2, 1:3, 1:4)))

  g2 <- graph.formula( s -+ a:b -+ t, a -+ 1:2:3 -+ b )
  cc <- stCuts(g2, source="s", target="t")
  expect_that(cc$cuts, equals(list(c(1,2), c(1,7), c(2,3,4,5,6),
                                   c(2,3,4,5,10), c(2,3,4,6,9),
                                   c(2,3,4,9,10), c(2,3,5,6,8),
                                   c(2,3,5,8,10), c(2,3,6,8,9),
                                   c(2,3,8,9,10), c(3,7))))
  expect_that(cc$partition1s,
              equals(list(1, c(1,3), c(1,2), c(1,2,7), c(1,2,6),
                          c(1,2,6,7), c(1,2,5), c(1,2,5,7), c(1,2,5,6),
                          c(1,2,5,6,7), c(1,2,5,6,7,3))))

  g3 <- graph.formula( s -+ a:b -+ t, a -+ 1:2:3:4:5 -+ b )
  cc <- stMincuts(g2, source="s", target="t")
  expect_that(cc$value, equals(2))
  expect_that(cc$cuts, equals(list(c(1,2), c(1,7), c(3,7))))
  expect_that(cc$partition1s, equals(list(1, c(1,3), c(1,3,2,7,6,5))))
})
