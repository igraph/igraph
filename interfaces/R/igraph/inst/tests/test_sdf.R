
context("sdf")

test_that("sdf works", {

  library(igraph)

  sdf <- igraph:::sdf(id=1:10, color="black")
  expect_that(as.data.frame(sdf),
              equals(data.frame(id=1:10, color="black")))

  ## access

  expect_that(sdf[1,"id"], equals(1))
  expect_that(sdf[1:4, "id"], equals(1:4))
  expect_that(sdf[, "id"], equals(1:10))

  expect_that(sdf[1, "color"], equals("black"))
  expect_that(sdf[1:4, "color"], equals(rep("black", 4)))
  expect_that(sdf[, "color"], equals(rep("black", 10)))

  ## set

  sdf2 <- sdf
  sdf2[5, "id"] <- 100
  expect_that(as.data.frame(sdf2),
              equals(data.frame(id=c(1:4,100,6:10), color="black")))

  sdf2 <- sdf
  sdf2[, "id"] <- 0
  expect_that(as.data.frame(sdf2),
              equals(data.frame(id=rep(0,10), color="black")))

  sdf2 <- sdf
  sdf2[2:10, "id"] <- 1
  expect_that(as.data.frame(sdf2),
              equals(data.frame(id=rep(1,10), color="black")))

  sdf2 <- sdf
  sdf2[, "color"] <- "white"
  expect_that(as.data.frame(sdf2),
              equals(data.frame(id=1:10, color="white")))

  sdf2 <- sdf
  sdf2[5:6, "color"] <- "white"
  expect_that(as.data.frame(sdf2),
              equals(data.frame(id=1:10, color=c(rep("black", 4),
                                           rep("white", 2),
                                           rep("black", 4)))))

})
