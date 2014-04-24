
context("scg")

## TODO: we only test that they run, not the results

test_that("SCG functions work", {

  library(igraph)

  tree <- graph.tree(10, 3, "undirected")
  treeM <- get.adjacency(tree, sparse=TRUE)
  treeM2 <- get.adjacency(tree, sparse=FALSE)

  args <- list(ev=1, nt=3, mtype="symmetric", algo="exact_scg",
               semproj=TRUE, epairs=TRUE)
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

  args[["ev"]] <- 3
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

  args[["ev"]] <- c(1,3)
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

###############################

  args <- list(ev=1, nt=2, mtype="stochastic", algo="exact_scg",
               semproj=TRUE, epairs=TRUE, stat.prob=TRUE)
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

  args[["ev"]] <- 3
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

  args[["ev"]] <- c(1,3)
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

###############################

  args <- list(ev=1, nt=2, mtype="laplacian", algo="exact_scg",
               semproj=TRUE, epairs=TRUE)
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

  args[["ev"]] <- 3
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

  args[["ev"]] <- c(1,3)
  do.call(scg, c(list(tree), args))
  do.call(scg, c(list(treeM), args))
  do.call(scg, c(list(treeM2), args))

})
