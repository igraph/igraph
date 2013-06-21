
context("walktrap.community")

test_that("walktrap.community works", {

  library(igraph)

  g <- graph.famous("Zachary")
  set.seed(42)
  wc <- walktrap.community(g)

  expect_that(modularity(g, membership(wc)), equals(modularity(wc)))
  expect_that(membership(wc), equals(c(1, 1, 2, 1, 5, 5, 5, 1, 2, 2,
                                       5, 1, 1, 2, 3, 3, 5, 1, 3, 1,
                                       3, 1, 3, 4, 4, 4, 3, 4, 2, 3,
                                       2, 2, 3, 3)))
  expect_that(length(wc), equals(5))
  expect_that(sizes(wc), equals(structure(c(9L, 7L, 9L, 4L, 5L), .Dim=5L,
    .Dimnames = structure(list(`Community sizes` = c("1", "2", "3", "4",
    "5")), .Names = "Community sizes"), class = "table")))

  d <- as.dendrogram(wc)
  expect_that(print(d), prints_text("2 branches.*34 members.*height 33"))
  expect_that(print(d[[1]]),
              prints_text("2 branches.*20 members.*height 31"))
  expect_that(print(d[[2]]),
              prints_text("2 branches.*14 members.*height 32"))
  m2 <- cutat(wc, no=3)
  expect_that(modularity(g, m2),
              equals(wc$modularity[length(wc$modularity)-2],
                     tolerance=1e-7))

})
