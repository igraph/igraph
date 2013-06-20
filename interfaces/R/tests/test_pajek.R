
context("Pajek file format")

test_that("writing Pajek files works", {

  library(igraph)

  g <- graph.ring(9)
  V(g)$color <- c("red", "green", "yellow")

  tc <- rawConnection(raw(0), "w")
  write.graph(g, format="pajek", file=tc)
  out <- rawToChar(rawConnectionValue(tc))
  close(tc)

  expect_that(out, equals("*Vertices 9\r\n1 \"1\" ic \"red\"\r\n2 \"2\" ic \"green\"\r\n3 \"3\" ic \"yellow\"\r\n4 \"4\" ic \"red\"\r\n5 \"5\" ic \"green\"\r\n6 \"6\" ic \"yellow\"\r\n7 \"7\" ic \"red\"\r\n8 \"8\" ic \"green\"\r\n9 \"9\" ic \"yellow\"\r\n*Edges\r\n1 2\r\n2 3\r\n3 4\r\n4 5\r\n5 6\r\n6 7\r\n7 8\r\n8 9\r\n1 9\r\n"))
  
})
