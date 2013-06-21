
context("is.chordal")

test_that("is.chordal works", {

  library(igraph)

  ## The examples from the Tarjan-Yannakakis paper
  g1 <- graph.formula(A-B:C:I, B-A:C:D, C-A:B:E:H, D-B:E:F,
                      E-C:D:F:H, F-D:E:G, G-F:H, H-C:E:G:I,
                      I-A:H)

  mc <- maximum.cardinality.search(g1)
  expect_that(mc, equals(list(alpha=c(9,4,6,8,3,5,7,2,1),
                              alpham1=c(9,8,5,2,6,3,7,4,1))))

  ic <- is.chordal(g1, fillin=TRUE)
  expect_that(ic$chordal, equals(FALSE))
  expect_that(unique(sort(ic$fillin)), equals(c(1,2,5,6,7,8)))
  expect_that(ic$newgraph, equals(NULL))

  g2 <- graph.formula(A-B:E, B-A:E:F:D, C-E:D:G, D-B:F:E:C:G,
                      E-A:B:C:D:F, F-B:D:E, G-C:D:H:I, H-G:I:J,
                      I-G:H:J, J-H:I)

  mc2 <- maximum.cardinality.search(g2)
  expect_that(mc2, equals(list(alpha=c(10,8,9,6,7,5,4,2,3,1),
                               alpham1=c(10,8,9,7,6,4,5,2,3,1))))
  
  ic2 <- is.chordal(g2, fillin=TRUE)
  expect_that(ic2, equals(list(chordal=TRUE, fillin=numeric(),
                               newgraph=NULL)))
})
