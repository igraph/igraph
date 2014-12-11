
context("Seeded graph matching")

test_that("SGM works", {
  library(igraph)
  set.seed(42)

  vc <- 10
  nos <- 3
  
  g1 <- erdos.renyi.game(vc, .5)
  randperm <- c(1:nos, nos + sample(vc-nos))
  g2 <- sample_correlated_gnp(g1, corr=.7, p=g1$p, perm=randperm)
  P  <-match_vertices (g1[], g2[], m=nos, start=matrix(1/(vc-nos), vc-nos, vc-nos),
            iteration=20)

  expect_that(c(1:nos, P$corr[,2]), equals(randperm))
  expect_that(apply(P$P != 0, 1, which), equals(randperm))
  expect_that(apply(P$D != 0, 1, which),
              equals(randperm[(nos+1):vc] - nos))

  ## Slightly bigger
  set.seed(42)

  vc <- 100
  nos <- 10

  g1 <- erdos.renyi.game(vc, .1);
  perm <- c(1:nos, sample(vc-nos)+nos)
  g2 <- sample_correlated_gnp(g1, corr=1, p=g1$p, perm=perm)

  P <- match_vertices(g1[], g2[], m=nos, start=matrix(1/(vc-nos), vc-nos, vc-nos),
           iteration=20)

  test_that(P$corr[,2], equals(perm[(nos+1):vc]))
  expect_that(apply(P$P != 0, 1, which), equals(perm))
  expect_that(apply(P$D != 0, 1, which),
              equals(perm[(nos+1):vc] - nos))
})

test_that("LSAP does not change input matrix", {

  x <- matrix(c(5, 1, 4, 3, 5, 2, 2, 4, 4), nrow = 3)
  solve_LSAP(x)

  expect_equal(x, matrix(c(5, 1, 4, 3, 5, 2, 2, 4, 4), nrow = 3))
  
})
