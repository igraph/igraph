
context("count.multiple")

test_that("count.multiple works", {
  library(igraph)
  set.seed(42)

  g <- barabasi.game(10, m=3, algorithm="bag")
  im <- is.multiple(g)
  cm <- count.multiple(g)
  expect_that(im, equals(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
                           FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
                           FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                           FALSE, FALSE, TRUE)))
  expect_that(cm, equals(c(3, 3, 3, 3, 3, 3, 1, 1, 1, 2, 1, 2, 1, 2,
                           2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2)))
  expect_that(count.multiple(simplify(g)),
              equals(rep(1, ecount(simplify(g)))))

  
  ## Direction of the edge is important
  expect_that(is.multiple(graph( c(1,2, 2,1) )), equals(c(FALSE, FALSE)))
  expect_that(is.multiple(graph( c(1,2, 2,1), dir=FALSE )),
              equals(c(FALSE, TRUE)))

  ## Remove multiple edges but keep multiplicity
  g <- barabasi.game(10, m=3, algorithm="bag")
  E(g)$weight <- 1
  g <- simplify(g)
  expect_that(any(is.multiple(g)), is_false())
  expect_that(E(g)$weight, equals(c(3, 2, 1, 2, 1, 3, 2, 1, 2, 1, 2,
                                    1, 1, 1, 1, 1, 1, 1)))

})
