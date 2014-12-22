
context("Detailed printing of vs and es")

test_that("vs printing", {

  set.seed(42)
  g <- make_graph(~ A - A:B:C, B - A:B:C) %>%
    set_vertex_attr("color", value = "red") %>%
    set_vertex_attr("weight", value = sample(1:10, 3))

  o1 <- c("+ 1/3 vertex, named:", "  name color weight",
          "1    A   red     10")
  expect_output(V(g)[[1]], paste(o1, collapse = "\n"), fixed = TRUE)

  o2 <- c("+ 1/3 vertex, named:", "  name color weight",
          "2    B   red      9")
  expect_output(V(g)[[2]], paste(o2, collapse = "\n"), fixed = TRUE)

  o3 <- c("+ 2/3 vertices, named:", "  name color weight",
          "1    A   red     10", "2    B   red      9")
  expect_output(V(g)[[1:2]], paste(o3, collapse = "\n"), fixed = TRUE)

  o4 <- c("+ 2/3 vertices, named:", "  name color weight",
          "2    B   red      9",  "3    C   red      3")
  expect_output(V(g)[[2:3]], paste(o4, collapse = "\n"), fixed = TRUE)

})

test_that("vs printing, complex attributes", {

  set.seed(42)
  g <- make_graph(~ A - A:B:C, B - A:B:C) %>%
    set_vertex_attr("color", value = "red") %>%
    set_vertex_attr("weight", value = sample(1:10, 3)) %>%
    set_vertex_attr("cplx", value = replicate(3, 1:4, simplify = FALSE))

  o1 <- c("+ 1/3 vertex, named:", "$name", "[1] \"A\"", "", "$color",
          "[1] \"red\"", "", "$weight", "[1] 10", "", "$cplx", "$cplx[[1]]",
          "[1] 1 2 3 4", "", "")
  expect_output(V(g)[[1]], paste(o1, collapse = "\n"), fixed = TRUE)

  o2 <- c("+ 2/3 vertices, named:", "$name", "[1] \"B\" \"C\"", "", "$color",
          "[1] \"red\" \"red\"", "", "$weight", "[1] 9 3", "", "$cplx",
          "$cplx[[1]]", "[1] 1 2 3 4", "", "$cplx[[2]]", "[1] 1 2 3 4",
          "", "")
  expect_output(V(g)[[2:3]], paste(o2, collapse = "\n"), fixed = TRUE)

})

test_that("es printing", {

  set.seed(42)
  g <- make_graph(~ A - A:B:C, B - A:B:C) %>%
    set_edge_attr("color", value = "red") %>%
    set_edge_attr("weight", value = sample(1:10, 3))

  o1 <- c("+ 1/3 edge (vertex names):",
          "  tail head tid hid color weight",
          "1    B    A   2   1   red     10")
  expect_output(E(g)[[1]], paste(o1, collapse = "\n"), fixed = TRUE)

  o2 <- c("+ 2/3 edges (vertex names):",
          "  tail head tid hid color weight",
          "2    C    A   3   1   red      9",
          "3    C    B   3   2   red      3")
  expect_output(E(g)[[2:3]], paste(o2, collapse = "\n"), fixed = TRUE)

})

test_that("es printing, complex attributes", {

  set.seed(42)
  g <- make_graph(~ A - A:B:C, B - A:B:C) %>%
    set_edge_attr("color", value = "red") %>%
    set_edge_attr("weight", value = sample(1:10, 3)) %>%
    set_edge_attr("cmpx", value = replicate(3, 1:4, simplify = FALSE))

  o1 <- c("+ 1/3 edge (vertex names):", "$color", "[1] \"red\"", "",
          "$weight", "[1] 10", "", "$cmpx", "$cmpx[[1]]", "[1] 1 2 3 4",
          "", "")
  expect_output(E(g)[[1]], paste(o1, collapse = "\n"), fixed = TRUE)

  o2 <- c("+ 2/3 edges (vertex names):", "$color", "[1] \"red\" \"red\"",
          "", "$weight", "[1] 9 3", "", "$cmpx", "$cmpx[[1]]",
          "[1] 1 2 3 4", "", "$cmpx[[2]]", "[1] 1 2 3 4", "", "")
  expect_output(E(g)[[2:3]], paste(o2, collapse = "\n"), fixed = TRUE)

})
