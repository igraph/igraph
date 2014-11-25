
context("VS/ES operators")

test_that("c on attached vs", {
  g <- make_ring(10)

  vg <- V(g)[1:5]
  vg2 <- V(g)[6:10]
  expect_equivalent(c(vg, vg2), V(g))
  expect_equal(get_vs_graph_id(c(vg, vg2)), get_graph_id(g))

  vg <- V(g)
  vg2 <- V(g)[FALSE]
  expect_equivalent(c(vg, vg2), V(g))
  expect_equivalent(c(vg2, vg), V(g))

  vg <- V(g)[c(2,5,6,8)]
  expect_equivalent(c(vg, vg), V(g)[c(2,5,6,8,2,5,6,8)])
})

test_that("c on detached vs", {
  g <- make_ring(10)

  vg <- V(g)[1:5]
  vg2 <- V(g)[6:10]

  vg3 <- V(g)
  vg4 <- V(g)[FALSE]

  vg5 <- V(g)[c(2,5,6,8)]
  vg6 <- V(g)[c(2,5,6,8,2,5,6,8)]

  rm(g)
  gc()

  expect_equivalent(c(vg, vg2), vg3)
  expect_equivalent(c(vg3, vg4), vg3)
  expect_equivalent(c(vg4, vg3), vg3)
  expect_equivalent(c(vg5, vg5), vg6)
})

test_that("c on attached vs, names", {
  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  vg <- V(g)[1:5]
  vg2 <- V(g)[6:10]
  expect_equivalent(c(vg, vg2), V(g))
  expect_equal(names(c(vg, vg2)), names(V(g)))

  vg <- V(g)
  vg2 <- V(g)[FALSE]
  expect_equivalent(c(vg, vg2), V(g))
  expect_equal(names(c(vg, vg2)), names(V(g)))
  expect_equivalent(c(vg2, vg), V(g))
  expect_equal(names(c(vg2, vg)), names(V(g)))

  vg <- V(g)[c(2,5,6,8)]
  expect_equivalent(c(vg, vg), V(g)[c(2,5,6,8,2,5,6,8)])
  expect_equal(names(c(vg, vg)), names(V(g)[c(2,5,6,8,2,5,6,8)]))
})

test_that("c on detached vs, names", {
  g <- make_ring(10)

  vg <- V(g)[1:5]
  vg2 <- V(g)[6:10]

  vg3 <- V(g)
  vg4 <- V(g)[FALSE]

  vg5 <- V(g)[c(2,5,6,8)]
  vg6 <- V(g)[c(2,5,6,8,2,5,6,8)]

  rm(g)
  gc()

  expect_equivalent(c(vg, vg2), vg3)
  expect_equal(names(c(vg, vg2)), names(vg3))
  expect_equivalent(c(vg3, vg4), vg3)
  expect_equal(names(c(vg3, vg4)), names(vg3))
  expect_equivalent(c(vg4, vg3), vg3)
  expect_equal(names(c(vg3, vg4)), names(vg3))
  expect_equivalent(c(vg5, vg5), vg6)
  expect_equal(names(c(vg5, vg5)), names(vg6))
})



test_that("union on attached vs", {

  g <- make_ring(10)

  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  vu <- union(v1, v2)
  expect_equivalent(vu, V(g))

  expect_equivalent(union(V(g)), V(g))

  v3 <- V(g)[FALSE]
  expect_equivalent(union(V(g), v3), V(g))
  expect_equivalent(union(v3, V(g), v3), V(g))
  expect_equivalent(union(v3), v3)
  expect_equivalent(union(v3, v3, v3), v3)
  expect_equivalent(union(v3, v3), v3)
})

test_that("union on detached vs", {

  g <- make_ring(10)

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  vu <- union(v1, v2)
  v3 <- V(g)[FALSE]

  rm(g)
  gc()

  expect_equivalent(vu, vg)

  expect_equivalent(union(vg), vg)

  expect_equivalent(union(vg, v3), vg)
  expect_equivalent(union(v3, vg, v3), vg)
  expect_equivalent(union(v3), v3)
  expect_equivalent(union(v3, v3, v3), v3)
  expect_equivalent(union(v3, v3), v3)

})

test_that("union on attached vs, names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  vu <- union(v1, v2)
  expect_equivalent(vu, V(g))
  expect_equal(names(vu), names(V(g)))

  expect_equivalent(union(V(g)), V(g))
  expect_equal(names(union(V(g))), names(V(g)))

  v3 <- V(g)[FALSE]
  expect_equivalent(union(V(g), v3), V(g))
  expect_equal(names(union(V(g), v3)), names(V(g)))

  expect_equivalent(union(v3, V(g), v3), V(g))
  expect_equal(names(union(v3, V(g), v3)), names(V(g)))

  expect_equivalent(union(v3), v3)
  expect_equal(names(union(v3)), names(v3))

  expect_equivalent(union(v3, v3, v3), v3)
  expect_equal(names(union(v3, v3, v3)), names(v3))

  expect_equivalent(union(v3, v3), v3)
  expect_equal(names(union(v3, v3)), names(v3))

})

test_that("union on detached vs, names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]

  rm(g)
  gc()

  vu <- union(v1, v2)
  expect_equivalent(vu, vg)
  expect_equal(names(vu), names(vg))

  expect_equivalent(union(vg), vg)
  expect_equal(names(union(vg)), names(vg))

  expect_equivalent(union(vg, v3), vg)
  expect_equal(names(union(vg, v3)), names(vg))

  expect_equivalent(union(v3, vg, v3), vg)
  expect_equal(names(union(v3, vg, v3)), names(vg))

  expect_equivalent(union(v3), v3)
  expect_equal(names(union(v3)), names(v3))

  expect_equivalent(union(v3, v3, v3), v3)
  expect_equal(names(union(v3, v3, v3)), names(v3))

  expect_equivalent(union(v3, v3), v3)
  expect_equal(names(union(v3, v3)), names(v3))

})



test_that("intersection on attached vs", {

  g <- make_ring(10)

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  v12 <- V(g)[6:7]
  v13 <- V(g)[FALSE]
  v14 <- V(g)[1:3]
  v24 <- V(g)[FALSE]

  vi1 <- intersection(v1, v2)
  expect_equivalent(vi1, v12)

  vi2 <- intersection(v1, v3)
  expect_equivalent(vi2, v13)

  vi3 <- intersection(v1, v4)
  expect_equivalent(vi3, v14)

  vi4 <- intersection(v1, vg)
  expect_equivalent(vi4, v1)

  vi5 <- intersection(v2, v4)
  expect_equivalent(vi5, v24)

  vi6 <- intersection(v3, vg)
  expect_equivalent(vi6, v3)

})

test_that("intersection on detached vs", {

  g <- make_ring(10)

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  v12 <- V(g)[6:7]
  v13 <- V(g)[FALSE]
  v14 <- V(g)[1:3]
  v24 <- V(g)[FALSE]

  rm(g)
  gc()

  vi1 <- intersection(v1, v2)
  expect_equivalent(vi1, v12)

  vi2 <- intersection(v1, v3)
  expect_equivalent(vi2, v13)

  vi3 <- intersection(v1, v4)
  expect_equivalent(vi3, v14)

  vi4 <- intersection(v1, vg)
  expect_equivalent(vi4, v1)

  vi5 <- intersection(v2, v4)
  expect_equivalent(vi5, v24)

  vi6 <- intersection(v3, vg)
  expect_equivalent(vi6, v3)

})

test_that("intersection on attached vs, names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  v12 <- V(g)[6:7]
  v13 <- V(g)[FALSE]
  v14 <- V(g)[1:3]
  v24 <- V(g)[FALSE]

  vi1 <- intersection(v1, v2)
  expect_equivalent(vi1, v12)
  expect_equal(names(vi1), names(v12))

  vi2 <- intersection(v1, v3)
  expect_equivalent(vi2, v13)
  expect_equal(names(vi2), names(v13))

  vi3 <- intersection(v1, v4)
  expect_equivalent(vi3, v14)
  expect_equal(names(vi3), names(v14))

  vi4 <- intersection(v1, vg)
  expect_equivalent(vi4, v1)
  expect_equal(names(vi4), names(v1))

  vi5 <- intersection(v2, v4)
  expect_equivalent(vi5, v24)
  expect_equal(names(vi5), names(v24))

  vi6 <- intersection(v3, vg)
  expect_equivalent(vi6, v3)
  expect_equal(names(vi6), names(v3))

})

test_that("intersection on detached vs, names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  v12 <- V(g)[6:7]
  v13 <- V(g)[FALSE]
  v14 <- V(g)[1:3]
  v24 <- V(g)[FALSE]

  rm(g)
  gc()

  vi1 <- intersection(v1, v2)
  expect_equivalent(vi1, v12)
  expect_equal(names(vi1), names(v12))

  vi2 <- intersection(v1, v3)
  expect_equivalent(vi2, v13)
  expect_equal(names(vi2), names(v13))

  vi3 <- intersection(v1, v4)
  expect_equivalent(vi3, v14)
  expect_equal(names(vi3), names(v14))

  vi4 <- intersection(v1, vg)
  expect_equivalent(vi4, v1)
  expect_equal(names(vi4), names(v1))

  vi5 <- intersection(v2, v4)
  expect_equivalent(vi5, v24)
  expect_equal(names(vi5), names(v24))

  vi6 <- intersection(v3, vg)
  expect_equivalent(vi6, v3)
  expect_equal(names(vi6), names(v3))

})



test_that("difference on attached vs", {

  g <- make_ring(10)

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  vr1 <- V(g)[8:10]
  vr2 <- V(g)
  vr3 <- V(g)[1:5]
  vr4 <- V(g)[4:7]
  vr5 <- V(g)[FALSE]
  vr6 <- V(g)[FALSE]

  vd1 <- difference(vg, v1)
  vd2 <- difference(vg, v3)
  vd3 <- difference(v1, v2)
  vd4 <- difference(v1, v4)
  vd5 <- difference(v3, v3)
  vd6 <- difference(v3, v4)

  expect_equivalent(vd1, vr1)
  expect_equivalent(vd2, vr2)
  expect_equivalent(vd3, vr3)
  expect_equivalent(vd4, vr4)
  expect_equivalent(vd5, vr5)
  expect_equivalent(vd6, vr6)

})

test_that("difference on detached vs", {

  g <- make_ring(10)

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  vr1 <- V(g)[8:10]
  vr2 <- V(g)
  vr3 <- V(g)[1:5]
  vr4 <- V(g)[4:7]
  vr5 <- V(g)[FALSE]
  vr6 <- V(g)[FALSE]

  rm(g)
  gc()

  vd1 <- difference(vg, v1)
  vd2 <- difference(vg, v3)
  vd3 <- difference(v1, v2)
  vd4 <- difference(v1, v4)
  vd5 <- difference(v3, v3)
  vd6 <- difference(v3, v4)

  expect_equivalent(vd1, vr1)
  expect_equivalent(vd2, vr2)
  expect_equivalent(vd3, vr3)
  expect_equivalent(vd4, vr4)
  expect_equivalent(vd5, vr5)
  expect_equivalent(vd6, vr6)

})

test_that("difference on attached vs, names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  vr1 <- V(g)[8:10]
  vr2 <- V(g)
  vr3 <- V(g)[1:5]
  vr4 <- V(g)[4:7]
  vr5 <- V(g)[FALSE]
  vr6 <- V(g)[FALSE]

  vd1 <- difference(vg, v1)
  vd2 <- difference(vg, v3)
  vd3 <- difference(v1, v2)
  vd4 <- difference(v1, v4)
  vd5 <- difference(v3, v3)
  vd6 <- difference(v3, v4)

  expect_equivalent(vd1, vr1)
  expect_equal(names(vd1), names(vr1))

  expect_equivalent(vd2, vr2)
  expect_equal(names(vd2), names(vr2))

  expect_equivalent(vd3, vr3)
  expect_equal(names(vd3), names(vr3))

  expect_equivalent(vd4, vr4)
  expect_equal(names(vd4), names(vr4))

  expect_equivalent(vd5, vr5)
  expect_equal(names(vd5), names(vr5))

  expect_equivalent(vd6, vr6)
  expect_equal(names(vd6), names(vr6))

})

test_that("difference on detached vs, names", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  vg <- V(g)
  v1 <- V(g)[1:7]
  v2 <- V(g)[6:10]
  v3 <- V(g)[FALSE]
  v4 <- V(g)[1:3]

  vr1 <- V(g)[8:10]
  vr2 <- V(g)
  vr3 <- V(g)[1:5]
  vr4 <- V(g)[4:7]
  vr5 <- V(g)[FALSE]
  vr6 <- V(g)[FALSE]

  rm(g)
  gc()

  vd1 <- difference(vg, v1)
  vd2 <- difference(vg, v3)
  vd3 <- difference(v1, v2)
  vd4 <- difference(v1, v4)
  vd5 <- difference(v3, v3)
  vd6 <- difference(v3, v4)

  expect_equivalent(vd1, vr1)
  expect_equal(names(vd1), names(vr1))

  expect_equivalent(vd2, vr2)
  expect_equal(names(vd2), names(vr2))

  expect_equivalent(vd3, vr3)
  expect_equal(names(vd3), names(vr3))

  expect_equivalent(vd4, vr4)
  expect_equal(names(vd4), names(vr4))

  expect_equivalent(vd5, vr5)
  expect_equal(names(vd5), names(vr5))

  expect_equivalent(vd6, vr6)
  expect_equal(names(vd6), names(vr6))

})



test_that("rev on attached vs", {

  for (i in 1:10) {
    g <- make_ring(10)
    idx <- seq_len(i)
    vg <- V(g)[idx]
    vgr <- V(g)[rev(idx)]
    vg2 <- rev(vg)
    expect_equivalent(vg2, vgr)
  }

})

test_that("rev on detached vs", {

  for (i in 1:10) {
    g <- make_ring(10)
    idx <- seq_len(i)
    vg <- V(g)[idx]
    vgr <- V(g)[rev(idx)]
    rm(g)
    gc()
    vg2 <- rev(vg)
    expect_equivalent(vg2, vgr)
  }

})

test_that("rev on attached vs, names", {

  for (i in 1:10) {
    g <- make_ring(10)
    V(g)$name <- letters[1:10]
    idx <- seq_len(i)
    vg <- V(g)[idx]
    vgr <- V(g)[rev(idx)]
    vg2 <- rev(vg)
    expect_equivalent(vg2, vgr)
    expect_equal(names(vg2), names(vgr))
  }

})

test_that("rev on detached vs, names", {

  for (i in 1:10) {
    g <- make_ring(10)
    V(g)$name <- letters[1:10]
    idx <- seq_len(i)
    vg <- V(g)[idx]
    vgr <- V(g)[rev(idx)]
    rm(g)
    gc()
    vg2 <- rev(vg)
    expect_equivalent(vg2, vgr)
    expect_equal(names(vg2), names(vgr))
  }

})

unique_tests <- list(
  list(1:5,        1:5),
  list(c(1,1,2:5), 1:5),
  list(c(1,1,1,1), 1),
  list(c(1,2,2,2), 1:2),
  list(c(2,2,1,1), 2:1),
  list(c(1,2,1,2), 1:2),
  list(c(),        c())
)

test_that("unique on attached vs", {

  sapply(unique_tests, function(d) {
    g <- make_ring(10)
    vg <- unique(V(g)[ d[[1]] ])
    vr <- V(g)[ d[[2]] ]
    expect_equivalent(vg, vr)
  })

})

test_that("unique on detached vs", {

  sapply(unique_tests, function(d) {
    g <- make_ring(10)
    vg <- V(g)[ d[[1]] ]
    vr <- V(g)[ d[[2]] ]
    rm(g)
    gc()
    vg <- unique(vg)
    expect_equivalent(vg, vr)
  })

})

test_that("unique on attached vs, names", {

  sapply(unique_tests, function(d) {
    g <- make_ring(10)
    V(g)$name <- letters[1:10]
    vg <- unique(V(g)[ d[[1]] ])
    vr <- V(g)[ d[[2]] ]
    expect_equivalent(vg, vr)
  })

})

test_that("unique on detached vs, names", {

  sapply(unique_tests, function(d) {
    g <- make_ring(10)
    V(g)$name <- letters[1:10]
    vg <- V(g)[ d[[1]] ]
    vr <- V(g)[ d[[2]] ]
    rm(g)
    gc()
    vg <- unique(vg)
    expect_equivalent(vg, vr)
  })

})
