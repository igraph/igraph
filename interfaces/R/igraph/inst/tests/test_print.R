
context("print.igraph")

test_that("print.igraph works", {

  library(igraph)
  igraph_options(print.full=TRUE)
  options(width=76)

  g <- make_ring(5)
  expect_that(summary(g), prints_text("attr:.* name[ ]*[(]g/c[)]"))
  expect_that(g, prints_text("attr:.* name[ ]*[(]g/c[)]"))
  expect_that(g, prints_text("1--2"))

  V(g)$name <- letters[1:vcount(g)]
  expect_that(summary(g), prints_text("name[ ]*[(]v/c[)]"))
  expect_that(g, prints_text("a--b"))

  set.seed(42)
  E(g)$weight <- sample(ecount(g))
  expect_that(summary(g), prints_text("weight[\n ]*[(]e/n[)]"))

  g$name <- "A ring"
  expect_that(summary(g), prints_text("A ring"))
  expect_that(print(g, v=T), prints_text("vertex attributes"))
  expect_that(print(g, e=T), prints_text("edges [(]vertex names[)] and"))

  set.seed(42)
  g2 <- sample_gnp(13, p=0.6, directed=TRUE)
  expect_that(g2, prints_text("1 ->"))

  g3 <- sample_gnp(20, p=0.8)
  expect_that(g3, prints_text("1 --"))

  g4 <- make_star(100)
  expect_that(g4, prints_text("2->1"))

  g5 <- make_star(100, mode="out")
  expect_that(g5, prints_text("1->"))

  g6 <- sample_pa(100, m=6, directed=FALSE)
  expect_that(g6, prints_text("     "))

  kite <- make_empty_graph(directed=FALSE) + LETTERS[1:10]
  kite <- kite + edges('A','B','A','C','A','D','A','F',
                       'B','D','B','E','B','G', 'C','D','C','F', 
                       'D','E','D','F','D','G', 'E','G', 
                       'F','G','F','H', 'G','H', 'H','I','I','J')
  expect_that(kite, prints_text("A -- "))
              
  igraph_options(print.full=FALSE)
})

test_that("print.igraph.es uses vertex names", {

  g <- make_directed_graph(c("A", "B"))
  expect_output(print(E(g)), "A\\s*->\\s*B")

})
