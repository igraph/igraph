
context("dominator_tree")

test_that("dominator_tree works", {
  library(igraph)
  g <- graph_from_literal(R-+A:B:C, A-+D, B-+A:D:E, C-+F:G, D-+L,
                 E-+H, F-+I, G-+I:J, H-+E:K, I-+K, J-+I,
                 K-+I:R, L-+H)
  dtree <- dominator_tree(g, root="R")

  dtree$dom <- V(g)$name[ as.vector(dtree$dom) ]
  dtree$leftout <- V(g)$name[ dtree$leftout ]
  expect_that(dtree$dom, equals(c("R", "R", "R", "R", "R", "C", "C",
                                  "D", "R", "R", "G", "R")))
  expect_that(dtree$leftout, equals(character()))
  expect_that(as_edgelist(dtree$domtree),
              equals(structure(c("R", "R", "R", "R", "R", "C", "C",
                                 "D", "R", "R", "G", "R", "A", "B",
                                 "C", "D", "E", "F", "G", "L", "H",
                                 "I", "J", "K"), .Dim = c(12L, 2L))))

})
