
context("Vertex and edge sequences")

test_that("we can create vertex/edge seqs", {

  g <- make_ring(10)
  V(g) %&&% expect_true(TRUE)
  E(g) %&&% expect_true(TRUE)

  V(g)$name <- letters[1:10]
  V(g) %&&% expect_true(TRUE)
  E(g) %&&% expect_true(TRUE)

  g <- make_ring(10)
  E(g)$name <- LETTERS[1:10]
  E(g) %&&% expect_true(TRUE)

})

test_that("vs/es keeps names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]
  vs <- V(g)

  expect_equal(vs$name, names(vs))

  vs2 <- vs[4:7]
  expect_equal(vs2$name, names(vs2))

  E(g)$name <- LETTERS[1:10]
  es <- E(g)

  expect_equal(es$name, names(es))

  es2 <- es[4:7]
  expect_equal(es2$name, names(es2))
})

test_that("vs/es refers to the graph", {

  g <- make_ring(10)
  vs <- V(g)
  es <- E(g)

  expect_identical(get_vs_graph(vs), g)
  expect_identical(get_es_graph(es), g)

})

test_that("vs/es refers to the original graph", {

  g <- g2 <- make_ring(10)
  vs <- V(g)
  es <- E(g)

  g <- g + 4

  expect_identical(get_vs_graph(vs), g2)
  expect_identical(get_es_graph(es), g2)

})

test_that("vs/es references are weak", {

  g <- make_ring(10)
  vs <- V(g)
  es <- E(g)

  rm(g)
  gc()

  expect_null(get_vs_graph(vs))
  expect_null(get_es_graph(es))
  
})

test_that("save/load breaks references", {

  g <- make_ring(10)
  vs <- V(g)
  es <- E(g)

  tmp <- tempfile()
  on.exit(try(unlink(tmp)))
  
  save(vs, es, file = tmp)
  rm(vs, es)
  gc()

  load(tmp)
  expect_null(get_vs_graph(vs))
  expect_null(get_es_graph(es))

})

test_that("we can use vs/es with broken refs", {

  g <- make_ring(10)
  vs <- V(g)
  es <- E(g)

  rm(g)
  gc()

  g2 <- make_ring(10)
  ## TODO
  
})

test_that("vs/es keeps names after graph is deleted", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]
  vs <- V(g)

  E(g)$name <- LETTERS[1:10]
  es <- E(g)
  
  rm(g)
  gc()
  
  expect_equal(names(vs), letters[1:10])

  vs2 <- vs[4:7]
  expect_equal(names(vs2), letters[4:7])

  expect_equal(names(es), LETTERS[1:10])

  es2 <- es[4:7]
  expect_equal(names(es2), LETTERS[4:7])
})

test_that("both edge and vertex names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]
  E(g)$name <- LETTERS[1:10]

  es <- E(g)
  expect_equal(as.vector(es), 1:10)
  expect_equal(names(es), LETTERS[1:10])
  el <- as_edgelist(g)
  expect_equal(attr(es, "vnames"), paste(el[,1], el[,2], sep = "|"))

  x1 <- es[LETTERS[4:7]]
  x2 <- E(g)[4:7]
  expect_equal(as.vector(x1), as.vector(x2))
  expect_equal(names(x1), names(x2))
  expect_equal(attr(x1, "vnames"), attr(x2, "vnames"))

  y1 <- es[c('a|b', 'd|e')]
  y2 <- E(g)[c(1,4)]
  expect_equal(as.vector(y1), as.vector(y2))
  expect_equal(names(y1), names(y2))
  expect_equal(attr(y1, "vnames"), attr(y2, "vnames"))

})

test_that("printing connected vs/es works", {

  g <- make_ring(10)
  vs <- V(g)
  es <- E(g)

  expect_output(print(vs), fixed = TRUE,
    "+ 10/10 vertices:\n [1]  1  2  3  4  5  6  7  8  9 10")
  expect_output(print(es), fixed = TRUE,
    "+ 10/10 edges:\n [1] 1-- 2 2-- 3 3-- 4 4-- 5 5-- 6 6-- 7 7-- 8 8-- 9 9--10 1--10")

  vs2 <- vs[1:5]
  es2 <- es[1:5]

  expect_output(print(vs2), fixed = TRUE,
    "+ 5/10 vertices:\n[1] 1 2 3 4 5")
  expect_output(print(es2), fixed = TRUE,
    "+ 5/10 edges:\n[1] 1--2 2--3 3--4 4--5 5--6")

  vs3 <- vs[numeric()]
  es3 <- es[numeric()]

  expect_output(print(vs3), fixed = TRUE, "+ 0/10 vertices:")
  expect_output(print(es3), fixed = TRUE, "+ 0/10 edges:")

  V(g)$name <- letters[1:10]
  vs <- V(g)
  es <- E(g)

  expect_output(print(vs), fixed = TRUE,
    "+ 10/10 vertices, named:\n [1] a b c d e f g h i j")
  expect_output(print(es), fixed = TRUE,
    "+ 10/10 edges (vertex names):\n [1] a--b b--c c--d d--e e--f f--g g--h h--i i--j a--j")

  vs2 <- vs[1:5]
  es2 <- es[1:5]

  expect_output(print(vs2), fixed = TRUE,
    "+ 5/10 vertices, named:\n[1] a b c d e")
  expect_output(print(es2), fixed = TRUE,
    "+ 5/10 edges (vertex names):\n[1] a--b b--c c--d d--e e--f")

  vs3 <- vs[numeric()]
  es3 <- es[numeric()]

  expect_output(print(vs3), fixed = TRUE, "+ 0/10 vertices, named:")
  expect_output(print(es3), fixed = TRUE, "+ 0/10 edges (vertex names):")
})

test_that("printing unconnected vs/es works", {

  g <- make_ring(10)
  vs <- V(g)
  es <- E(g)

  rm(g)
  gc()

  expect_output(print(vs), fixed = TRUE,
    "+ 10/? vertices:\n [1]  1  2  3  4  5  6  7  8  9 10")                
  expect_output(print(es), fixed = TRUE,
    "+ 10/? edges, unknown graph:\n [1]  1  2  3  4  5  6  7  8  9 10")

  g <- make_ring(10)
  V(g)$name <- letters[1:10]
  vs <- V(g)
  es <- E(g)

  rm(g)
  gc()

  expect_output(print(vs), fixed = TRUE,
    "+ 10/? vertices, named:\n [1] a b c d e f g h i j")
  expect_output(print(es), fixed = TRUE,
    "+ 10/? edges, unknown graph (vertex names):\n [1] a|b b|c c|d d|e e|f f|g g|h h|i i|j a|j")
  
})
