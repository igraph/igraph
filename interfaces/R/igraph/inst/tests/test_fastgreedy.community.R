
context("fastgreedy.community")

test_that("fastgreedy.community works", {

  library(igraph)
  set.seed(42)

  g <- graph.famous("Zachary")
  fc <- fastgreedy.community(g)

  expect_that(modularity(g, fc$membership), equals(max(fc$modularity)))
  expect_that(membership(fc), equals(c(1, 3, 3, 3, 1, 1, 1, 3, 2, 3,
                                       1, 1, 3, 3, 2, 2, 1, 3, 2, 1,
                                       2, 3, 2, 2, 2, 2, 2, 2, 2, 2,
                                       2, 2, 2, 2)))
  expect_that(length(fc), equals(3))
  expect_that(as.numeric(sizes(fc)), equals(c(8, 17, 9)))

  d <- as.dendrogram(fc)
  expect_that(print(d), prints_text("2 branches.*34 members.*height 33"))
  expect_that(print(d[[1]]),
              prints_text("2 branches.*17 members.*height 32"))
  expect_that(print(d[[2]]),
              prints_text("2 branches.*17 members.*height 30"))
  m2 <- cutat(fc, no=3)
  expect_that(modularity(g, m2),
              equals(fc$modularity[length(fc$modularity)-2]))
})
