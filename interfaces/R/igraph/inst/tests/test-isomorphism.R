
context("New isomorphism API")

test_that("isomorphic", {

  g <- graph_(from_literal(A - B - C - A))
  expect_true(isomorphic(g, g))
  expect_true(isomorphic(g, g, method = "direct"))
  expect_true(isomorphic(g, g, method = "vf2"))
  expect_true(isomorphic(g, g, method = "bliss"))

  g2 <- graph_(from_literal(A - B - C))
  expect_false(isomorphic(g, g2))
  expect_false(isomorphic(g, g2, method = "direct"))
  expect_false(isomorphic(g, g2, method = "vf2"))
  expect_false(isomorphic(g, g2, method = "bliss"))
  
})

test_that("subgraph_isomorphic", {

  g <- graph_(from_literal(A - B - C - D - E - A))
  g2 <- graph_(from_literal(A - B - C - D))

  expect_true(subgraph_isomorphic(g2, g))
  expect_true(subgraph_isomorphic(g2, g, method = "vf2"))
  expect_true(subgraph_isomorphic(g2, g, method = "lad"))

  g3 <- graph_(from_literal(A - B - C - A))
  expect_false(subgraph_isomorphic(g3, g))
  expect_false(subgraph_isomorphic(g3, g, method = "vf2"))
  expect_false(subgraph_isomorphic(g3, g, method = "lad"))
  
})

test_that("count_isomorphisms", {

  g <- graph_(from_literal(A - B - C - D - A))
  expect_equal(count_isomorphisms(g, g), 8)

  g2 <- graph_(from_literal(A - B - C - A))
  expect_equal(count_isomorphisms(g, g2), 0)
  
})

test_that("count_subgraph_isomorphisms", {

  g <- graph_(from_literal(A - B - C - D - A))
  g2 <- graph_(from_literal(A - B - C - D))

  expect_equal(count_subgraph_isomorphisms(g2, g, method = "lad"), 8)
  expect_equal(count_subgraph_isomorphisms(g2, g, method = "vf2"), 8)

  g3 <- graph_(from_literal(A - B - C - A))
  expect_equal(count_subgraph_isomorphisms(g3, g, method = "lad"), 0)
  expect_equal(count_subgraph_isomorphisms(g3, g, method = "vf2"), 0)
  
})

test_that("isomorphisms", {

  g <- graph_(from_literal(A - B - C - D - A))
  g2 <- graph_(from_literal(W - X - Y - Z - W))

  res <- list(V(g2)[1,2,3,4],
              V(g2)[1,4,3,2],
              V(g2)[2,1,4,3],
              V(g2)[2,3,4,1],
              V(g2)[3,2,1,4],
              V(g2)[3,4,1,2],
              V(g2)[4,1,2,3],
              V(g2)[4,3,2,1])
                
  expect_equivalent(isomorphisms(g, g2), res)

  g3 <- graph_(from_literal(X - Y - Z - X))
  expect_equal(isomorphisms(g, g3), list())

})

test_that("subgraph_isomorphisms, lad", {

  g <- graph_(from_literal(A - B - C - D - A))
  g2 <- graph_(from_literal(Z - X - Y))

  res <- list(V(g)[1,4,3],
              V(g)[1,2,3],
              V(g)[2,1,4],
              V(g)[2,3,4],
              V(g)[3,2,1],
              V(g)[3,4,1],
              V(g)[4,3,2],
              V(g)[4,1,2])
  
  expect_equivalent(subgraph_isomorphisms(g2, g, method = "lad"), res)

  g3 <- graph_(from_literal(X - Y - Z - X))
  expect_equal(subgraph_isomorphisms(g3, g, method = "lad"), list())

})

test_that("subgraph_isomorphisms, vf2", {

  g <- graph_(from_literal(A - B - C - D - A))
  g2 <- graph_(from_literal(Z - X - Y))

  res <- list(V(g)[1,2,3],
              V(g)[1,4,3],
              V(g)[2,1,4],
              V(g)[2,3,4],
              V(g)[3,2,1],
              V(g)[3,4,1],
              V(g)[4,1,2],
              V(g)[4,3,2])
  
  expect_equivalent(subgraph_isomorphisms(g2, g, method = "vf2"), res)

  g3 <- graph_(from_literal(X - Y - Z - X))
  expect_equal(subgraph_isomorphisms(g3, g, method = "vf2"), list())

})
