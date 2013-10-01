
context("motifs")

test_that("motif finding works", {

  library(igraph)
  set.seed(123)

  b <- erdos.renyi.game(10000, 4/10000, directed=TRUE)

  mno <- graph.motifs.no(b)

  mno0 <- graph.motifs.no(b, cut.prob=c(1/3, 0, 0))
  mno1 <- graph.motifs.no(b, cut.prob=c(0, 0, 1/3))
  mno2 <- graph.motifs.no(b, cut.prob=c(0, 1/3, 0))
  expect_that(c(mno0/mno, mno1/mno, mno2/mno),
              equals(c(0.654821903845065, 0.666289144345659,
                       0.668393831285275)))

  mno3 <- graph.motifs.no(b, cut.prob=c(0, 1/3, 1/3))
  mno4 <- graph.motifs.no(b, cut.prob=c(1/3, 0, 1/3))
  mno5 <- graph.motifs.no(b, cut.prob=c(1/3, 1/3, 0))
  expect_that(c(mno3/mno, mno4/mno, mno5/mno),
              equals(c(0.443959957465819, 0.441952797125797,
                       0.446004870037941) ))

######################

  set.seed(123)
  b <- erdos.renyi.game(10000, 4/10000, directed=TRUE)

  m <- graph.motifs(b)

  m0 <- graph.motifs(b, cut.prob=c(1/3, 0, 0))
  m1 <- graph.motifs(b, cut.prob=c(0, 1/3, 0))
  m2 <- graph.motifs(b, cut.prob=c(0, 0, 1/3))
  expect_that(m0/m, equals(c(NA, NA, 0.653972107372707, NA,
                             0.653993015279859, 0.612244897959184,
                             0.657514670174019, 0.63013698630137, NaN,
                             0.538461538461538, NaN,
                             0.565217391304348, NaN, NaN, NaN, NaN)))
  expect_that(m1/m, equals(c(NA, NA, 0.669562138856225, NA,
                             0.66808158454082, 0.73469387755102,
                             0.670819000404694, 0.657534246575342,
                             NaN, 0.769230769230769, NaN,
                             0.739130434782609, NaN, NaN, NaN, NaN) ))
  expect_that(m2/m, equals(c(NA, NA, 0.666451718949538, NA,
                             0.665291458452201, 0.591836734693878,
                             0.666683528935654, 0.671232876712329,
                             NaN, 0.753846153846154, NaN,
                             0.565217391304348, NaN, NaN, NaN, NaN) ))

  m3 <- graph.motifs(b, cut.prob=c(0, 1/3, 1/3))
  m4 <- graph.motifs(b, cut.prob=c(1/3, 1/3, 0))
  m5 <- graph.motifs(b, cut.prob=c(1/3, 1/3, 0))
  expect_that(m3/m, equals(c(NA, NA, 0.445611905574732, NA,
                             0.442789875290769, 0.448979591836735,
                             0.444695973290166, 0.424657534246575,
                             NaN, 0.369230769230769, NaN,
                             0.608695652173913, NaN, NaN, NaN, NaN)))

  expect_that(m4/m, equals(c(NA, NA, 0.439251981944392, NA,
                             0.439284975327761, 0.73469387755102,
                             0.445088021044112, 0.465753424657534,
                             NaN, 0.630769230769231, NaN,
                             0.565217391304348, NaN, NaN, NaN, NaN) ))

  expect_that(m5/m, equals(c(NA, NA, 0.439985332979302, NA,
                             0.440288166730411, 0.346938775510204,
                             0.44159753136382, 0.452054794520548, NaN,
                             0.323076923076923, NaN,
                             0.347826086956522, NaN, NaN, NaN, NaN) ))
})
