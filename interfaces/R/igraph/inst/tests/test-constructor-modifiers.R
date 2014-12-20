
context("Constructor modifiers")


test_that("without_attr", {

  set.seed(42)
  g <- sample_gnp(10, 2/10) %>%
    delete_graph_attr("name") %>%
    delete_graph_attr("type") %>%
    delete_graph_attr("loops") %>%
    delete_graph_attr("p")
  
  set.seed(42)
  g2 <- sample_(gnp(10, 2/10), without_attr())

  expect_equivalent(g, g2)
  expect_equal(graph_attr_names(g2), character())
  expect_equal(graph_vertex_names(g2), character())
  expect_equal(graph_edge_names(g2), character())

})


test_that("without_loops", {

  g <- make_graph(~ A - A:B:C, B - A:B:C, simplify = FALSE) %>%
     simplify(remove.multiple = FALSE)

  g2 <- make_(from_literal(A - A:B:C, B - A:B:C, simplify = FALSE),
              without_loops())
              
  expect_equivalent(g, g2)
  expect_true(all(!which_loop(g2)))

})


test_that("without_multiple", {

  g <- make_graph(~ A - A:B:C, B - A:B:C, simplify = FALSE) %>%
     simplify(remove.loops = FALSE)

  g2 <- make_(from_literal(A - A:B:C, B - A:B:C, simplify = FALSE),
              without_multiples())
              
  expect_equivalent(g, g2)
  expect_true(all(!which_multiple(g2)))

})


test_that("simplified", {

  g <- make_graph(~ A - A:B:C, B - A:B:C)

  g2 <- make_(from_literal(A - A:B:C, B - A:B:C, simplify = FALSE),
              simplified())
  
  expect_equivalent(g, g2)
  expect_true(all(!which_multiple(g2)))
  expect_true(all(!which_loop(g2)))
  
})


test_that("with_vertex_", {

  g <- make_graph(~ A - A:B:C, B - A:B:C) %>%
    set_vertex_attr("color", value = "red") %>%
    set_vertex_attr("foo", value = paste0("xx", 1:3))

  g2 <- make_(from_literal(A - A:B:C, B - A:B:C),
              with_vertex_(color = "red",
                           foo = paste0("xx", 1:3))
              )

  expect_equivalent(g, g2)
  expect_equal(V(g2)$color, rep("red", gorder(g2)))
  expect_equal(V(g2)$foo, paste0("xx", 1:3))
  
})


test_that("with_edge_", {

  g <- make_graph(~ A - A:B:C, B - A:B:C) %>%
    set_edge_attr("color", value = "red") %>%
    set_edge_attr("foo", value = seq_len(3))

  g2 <- make_(from_literal(A - A:B:C, B - A:B:C),
              with_edge_(color = "red",
                         foo = seq_len(3)))

  expect_equivalent(g, g2)
  expect_equal(E(g)$color, E(g2)$color)
  expect_equal(E(g)$foo, E(g2)$foo)
  
})


test_that("with_graph_", {

  g <- make_graph(~ A - A:B:C, B - A:B:C) %>%
    set_graph_attr("color", value = "red") %>%
    set_graph_attr("foo", value = 1:5)

  g2 <- make_(from_literal(A - A:B:C, B - A:B:C),
              with_graph_(color = "red",
                          foo = 1:5))

  expect_equivalent(g, g2)
  expect_equal(g$color, g2$color)
  expect_equal(g$foo, g2$foo)
  
})
