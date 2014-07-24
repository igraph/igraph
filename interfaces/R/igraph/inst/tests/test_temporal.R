
context("Temporal graph")

test_that("We can create temporal graphs", {

  library(igraph)
  g <- graph.temporal(1:10, v_on = 1, e_on = 1:5, calendar = 1:20)

  expect_that(is.temporal(g), is_true())

  expect_that(g$calendar, equals(1:20))
  expect_that(g$now, equals(Inf))

  g$now <- 5
  expect_that(g$now, equals(5))

  g$now <- Inf
  expect_that(g$now, equals(Inf))

  expect_that(g$now <- 21, throws_error("Unknown time point"))

  #####

  g2 <- graph.temporal(1:10, v_on = 1,
                      e_on = c("2001", "2001", "2010", "2010", "2014"))
  g2$now <- "2010"
  expect_that(g2$now, equals("2010"))

  expect_that(g2$now <- "2005", throws_error("Unknown time point"))

  expect_that(is.temporal(g2), is_true())

  #####

  g3 <- graph(1:10)
  g3 <- as.temporal(g3)
  expect_that(is.temporal(g3), is_true())
  expect_that(g3$calendar, equals(1))

})

test_that("on() and off() works", {

  library(igraph)
  g <- graph.temporal(1:10, v_on = 1, e_on = 1:5, calendar = 1:20)

  expect_that(on(V(g)), equals(rep(1, 10)))
  expect_that(on(E(g)), equals(1:5))
  expect_that(off(V(g)), equals(rep(Inf, 10)))
  expect_that(off(E(g)), equals(rep(Inf, 5)))

  vr <- vertex_range(g)
  er <- edge_range(g)

  expect_that(vr, equals(list(on=rep(1, 10), off=rep(Inf, 10))))
  expect_that(er, equals(list(on=1:5, off=rep(Inf, 5))))

  g$now <- 2

  expect_that(on(V(g)), equals(rep(1, 10)))
  expect_that(on(E(g)), equals(1:2))
  expect_that(off(V(g)), equals(rep(Inf, 10)))
  expect_that(off(E(g)), equals(rep(Inf, 2)))

  vr <- vertex_range(g)
  er <- edge_range(g)

  expect_that(vr, equals(list(on=rep(1, 10), off=rep(Inf, 10))))
  expect_that(er, equals(list(on=1:5, off=rep(Inf, 5))))

  ## Symbolic calendar

  e_on <- c("2001", "2001", "2005", "2010", "2013")
  g <- graph.temporal(1:10, v_on = "2001", e_on = e_on)

  expect_that(on(V(g)), equals(rep("2001", 10)))
  expect_that(on(E(g)), equals(e_on))
  expect_that(off(V(g)), equals(rep(Inf, 10)))
  expect_that(off(E(g)), equals(rep(Inf, 5)))

  vr <- vertex_range(g)
  er <- edge_range(g)

  expect_that(vr, equals(list(on=rep("2001", 10), off=rep(Inf, 10))))
  expect_that(er, equals(list(on=e_on, off=rep(Inf, 5))))

  g$now <- "2005"

  expect_that(on(V(g)), equals(rep("2001", 10)))
  expect_that(on(E(g)), equals(e_on[1:3]))
  expect_that(off(V(g)), equals(rep(Inf, 10)))
  expect_that(off(E(g)), equals(rep(Inf, 3)))

  vr <- vertex_range(g)
  er <- edge_range(g)

  expect_that(vr, equals(list(on=rep("2001", 10), off=rep(Inf,10))))
  expect_that(er, equals(list(on=e_on, off=rep(Inf, 5))))

})
