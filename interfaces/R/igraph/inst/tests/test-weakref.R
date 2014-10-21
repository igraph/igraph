
context("Weak references")

test_that("we can create weak references", {

  g <- new.env()
  g$foo <- "bar"
  value <- "foobar"
  vs <- make_weak_ref(key = g, value = value)

  expect_identical(typeof(vs), "weakref")
  expect_identical(weak_ref_key(vs), g)
  expect_identical(weak_ref_value(vs), value)

})

test_that("weak references are weak", {

  g <- new.env()
  g$foo <- "bar"
  value <- "foobar"
  vs <- make_weak_ref(key = g, value = value)

  rm(g)
  gc()
  expect_null(weak_ref_key(vs))
  expect_null(weak_ref_value(vs))

})

test_that("weak reference finalizer is called", {

  g <- new.env()
  g$foo <- "bar"
  value <- "foobar"
  hello <- ""
  fin <- function(env) hello <<- "world"
  vs <- make_weak_ref(key = g, value = value, finalizer = fin)

  rm(g)
  gc()

  expect_equal(hello, "world")

})

test_that("weak reference on an embedded env", {

  g <- list(yes = new.env())
  g[[1]]$foo <- "bar"
  value <- "foobar"
  vs <- make_weak_ref(key = g[[1]], value = value)

  rm(g)
  gc()
  expect_null(weak_ref_key(vs))
  expect_null(weak_ref_value(vs))
})

test_that("embed myself, and weak ref", {

  g <- list(yes = new.env())
  assign("foo", g, envir = g[[1]])
  value <- "foobar"
  hello <- ""
  fin <- function(env) hello <<- "world"
  vs <- make_weak_ref(key = g[[1]], value = value, finalizer = fin)

  rm(g)
  gc()
  expect_null(weak_ref_key(vs))
  expect_null(weak_ref_value(vs))
  expect_equal(hello, "world")

})

test_that("embed myself, and weak ref as attribute", {

  g <- list(yes = new.env())
  assign("foo", g, envir = g[[1]])
  value <- "foobar"
  hello <- ""
  fin <- function(env) hello <<- "world"
  z <- "footoo"
  attr(z, "env") <- make_weak_ref(key = g[[1]], value = value,
                                  finalizer = fin)

  rm(g)
  gc()
  expect_null(weak_ref_key(attr(z, "env")))
  expect_null(weak_ref_value(attr(z, "env")))
  expect_equal(hello, "world")

})

test_that("weak refs work for vs", {

  g <- make_ring(10)
  vs <- V(g)
  expect_true(!is.null(get_vs_ref(g)))
  expect_true(!is.null(weak_ref_key(attr(vs, "env"))))

  rm(g)
  gc()
  expect_null(weak_ref_key(attr(vs, "env")))
})
