
context("Random walks")

test_that("Undirected", {

  set.seed(42)
  g <- make_ring(10)
  w <- random_walk(g, start = 1, steps = 10)
  expect_equivalent(w, structure(c(1L, 10L, 9L, 8L, 9L, 10L, 9L, 10L,
                                   1L, 10L), class = "igraph.vs"))

})

test_that("Directed", {

  set.seed(42)
  g <- make_ring(10, directed = TRUE)
  w <- as_ids(random_walk(g, start = 1, steps = 5))
  expect_equal(w, 1:5)

  w2 <- as_ids(random_walk(g, start = 5, steps = 5, mode = "in"))
  expect_equal(w2, 5:1)

  set.seed(42)
  w3 <- as_ids(random_walk(g, start = 1, steps = 5, mode = "all"))
  expect_equal(w3, c(1, 10, 9, 8, 9))

})
