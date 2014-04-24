
context("graph.data.frame")

test_that("graph.data.frame works", {

  library(igraph) ; igraph.options(print.full=TRUE)

  actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                         "Esmeralda"),
                       age=c(48,33,45,34,21),
                       gender=c("F","M","F","M","F"),
                       stringsAsFactors=FALSE)
  relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                            "David", "Esmeralda"),
                          to=c("Alice", "Bob", "Alice", "Alice",
                            "Bob", "Alice"),
                          same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                          friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3),
                          stringsAsFactors=FALSE)
  g <- graph.data.frame(relations, directed=TRUE, vertices=actors)

  df <- get.data.frame(g, what="both")
  expect_that(df$vertices, is_equivalent_to(actors))
  expect_that(df$edges, equals(relations))

})

test_that("graph.data.frame works on matrices", {

  library(igraph)

  el <- cbind(1:5,5:1,weight=1:5)
  g <- graph.data.frame(el)
  g <- remove.vertex.attribute(g, "name")
  el2 <- get.data.frame(g)
  expect_that(as.data.frame(el), is_equivalent_to(el2))

})
