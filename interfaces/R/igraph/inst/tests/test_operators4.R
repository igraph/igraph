
context("operators on named graphs")


test_that("disjoint union works for named graphs", {

  library(igraph)

  g1 <- g2 <- graph.ring(10)
  g1$foo <- "bar"
  V(g1)$name <- letters[ 1:10]
  V(g2)$name <- letters[11:20]
  E(g1)$weight <- 1:10
  E(g2)$weight <- 10:1

  V(g1)$a1 <- 1:10
  V(g2)$a2 <- 11:20
  
  E(g1)$b1 <- 1:10
  E(g2)$b2 <- 11:20
  
  g <- g1 + g2

  expect_that(sort(list.graph.attributes(g)),
              equals(c("circular_1", "circular_2", "foo", "mutual_1",
                       "mutual_2", "name_1", "name_2")))
  expect_that(sort(list.vertex.attributes(g)),
              equals(c("a1", "a2", "name")))
  expect_that(sort(list.edge.attributes(g)),
              equals(c("b1", "b2", "weight")))

  expect_that(V(g)$name, equals(letters[1:20]))
  expect_that(V(g)$a1, equals(c(1:10, rep(NA, 10))))
  expect_that(V(g)$a2, equals(c(rep(NA, 10), 11:20)))

  expect_that(E(g)$weight, equals(c(1:10, 10:1)))
  expect_that(E(g)$b1, equals(c(1:10, rep(NA, 10))))
  expect_that(E(g)$b2, equals(c(rep(NA, 10), 11:20)))

})

test_that("union of unnamed graphs works", {

  library(igraph)

  g1 <- graph.ring(10)
  g2 <- graph.ring(13)
  g1$foo <- "bar"
  E(g1)$weight <- 1:10
  E(g2)$weight <- 13:1
  
  V(g1)$a1 <- 1:10
  V(g2)$a2 <- 11:23
  
  E(g1)$b1 <- letters[1:10]
  E(g2)$b2 <- letters[11:23]
  
  g <- graph.union(g1, g2)

  expect_that(sort(list.graph.attributes(g)),
              equals(c("circular_1", "circular_2", "foo", "mutual_1",
                       "mutual_2", "name_1", "name_2")))
  expect_that(sort(list.vertex.attributes(g)),
              equals(c("a1", "a2")))
  expect_that(sort(list.edge.attributes(g)),
              equals(c("b1", "b2", "weight_1", "weight_2")))

  df1 <- get.data.frame(g)
  df1 <- df1[ order(df1$from, df1$to), c(1,2,3,5,4,6)]
  df2 <- merge(get.data.frame(g1), get.data.frame(g2),
               by=c("from", "to"), all=TRUE)
  rownames(df1) <- seq_len(nrow(df1))
  colnames(df2) <- c("from", "to", "weight_1", "b1", "weight_2", "b2")
  expect_that(df1, equals(df2))

})

test_that("union of named graphs works", {

  library(igraph)

  g1 <- graph.ring(10)
  g2 <- graph.ring(13)
  V(g1)$name <- letters[seq_len(vcount(g1))]
  V(g2)$name <- letters[seq_len(vcount(g2))]

  g1$foo <- "bar"
  E(g1)$weight <- 1:10
  E(g2)$weight <- 13:1
  
  V(g1)$a1 <- 1:10
  V(g2)$a2 <- 11:23

  E(g1)$b1 <- letters[1:10]
  E(g2)$b2 <- letters[11:23]

  g <- graph.union(g1, g2)

  expect_that(sort(list.graph.attributes(g)),
              equals(c("circular_1", "circular_2", "foo",
                       "mutual_1", "mutual_2", "name_1", "name_2")))
  expect_that(sort(list.vertex.attributes(g)),
              equals(c("a1", "a2", "name")))
  expect_that(sort(list.edge.attributes(g)),
              equals(c("b1", "b2", "weight_1", "weight_2")))
  
  df1 <- get.data.frame(g, what="both")

  g.v <- read.table(stringsAsFactors=FALSE, textConnection("
  a1 a2 name
a  1 11    a
b  2 12    b
c  3 13    c
d  4 14    d
e  5 15    e
f  6 16    f
g  7 17    g
h  8 18    h
i  9 19    i
j 10 20    j
k NA 21    k
l NA 22    l
m NA 23    m
"))
  expect_that(df1$vertices, equals(g.v))

  g.e <- read.table(stringsAsFactors=FALSE, textConnection("
   from to weight_1 weight_2   b1   b2
1     l  m       NA        2   NA    v
2     k  l       NA        3   NA    u
3     j  k       NA        4   NA    t
4     i  j        9        5    i    s
5     h  i        8        6    h    r
6     g  h        7        7    g    q
7     f  g        6        8    f    p
8     e  f        5        9    e    o
9     d  e        4       10    d    n
10    c  d        3       11    c    m
11    b  c        2       12    b    l
12    a  m       NA        1   NA    w
13    a  j       10       NA    j   NA
14    a  b        1       13    a    k
"))
  rownames(df1$edges) <- rownames(df1$edges)
  expect_that(df1$edges, equals(g.e))
  
})

test_that("intersection of named graphs works", {

  library(igraph)

  g1 <- graph.ring(10)
  g2 <- graph.ring(13)
  V(g1)$name <- letters[V(g1)]
  V(g2)$name <- letters[V(g2)]

  g1$foo <- "bar"
  E(g1)$weight <- 1:10
  E(g2)$weight <- 13:1

  V(g1)$a1 <- 1:10
  V(g2)$a2 <- 11:23

  E(g1)$b1 <- letters[1:10]
  E(g2)$b2 <- letters[11:23]

  g <- graph.intersection(g1, g2)

  expect_that(sort(list.graph.attributes(g)),
              equals(c("circular_1", "circular_2", "foo", "mutual_1",
                       "mutual_2", "name_1", "name_2")))
  expect_that(sort(list.vertex.attributes(g)),
              equals(c("a1", "a2", "name")))
  expect_that(sort(list.edge.attributes(g)),
              equals(c("b1", "b2", "weight_1", "weight_2")))

  df1 <- get.data.frame(g, what="both")

  g.e <- read.table(stringsAsFactors=FALSE, textConnection("
  from to weight_1 weight_2 b1 b2
1    i  j        9        5  i  s
2    h  i        8        6  h  r
3    g  h        7        7  g  q
4    f  g        6        8  f  p
5    e  f        5        9  e  o
6    d  e        4       10  d  n
7    c  d        3       11  c  m
8    b  c        2       12  b  l
9    a  b        1       13  a  k
"))
  rownames(df1$edges) <- rownames(df1$edges)
  expect_that(df1$edges, equals(g.e))

  g.v <- read.table(stringsAsFactors=FALSE, textConnection("
  a1 a2 name
a  1 11    a
b  2 12    b
c  3 13    c
d  4 14    d
e  5 15    e
f  6 16    f
g  7 17    g
h  8 18    h
i  9 19    i
j 10 20    j
"))
  expect_that(df1$vertices, equals(g.v))

  gg <- graph.intersection(g1, g2, keep.all.vertices=TRUE)

  df2 <- get.data.frame(gg, what="both")

  rownames(df2$edges) <- rownames(df2$edges)
  expect_that(df2$edges, equals(g.e))

  gg.v <- read.table(stringsAsFactors=FALSE, textConnection("
  a1 a2 name
a  1 11    a
b  2 12    b
c  3 13    c
d  4 14    d
e  5 15    e
f  6 16    f
g  7 17    g
h  8 18    h
i  9 19    i
j 10 20    j
k NA 21    k
l NA 22    l
m NA 23    m
"))
  expect_that(df2$vertices, equals(gg.v))
})

test_that("difference of named graphs works", {

  library(igraph)

  g1 <- graph.ring(10)
  g2 <- graph.star(11, center=11, mode="undirected")
  V(g1)$name <- letters[1:10]
  V(g2)$name <- letters[1:11]
  g <- g1 %u% g2

  sg <- graph.ring(4)
  V(sg)$name <- letters[c(1,2,3,11)]

  df1 <- get.data.frame(g - sg, what="both")

  t1.e <- read.table(stringsAsFactors=FALSE,
                                           textConnection("
   from to
1     a  j
2     b  k
3     c  d
4     j  k
5     i  k
6     h  k
7     g  k
8     f  k
9     e  k
10    d  k
11    d  e
12    e  f
13    f  g
14    g  h
15    h  i
16    i  j
"))
  rownames(df1$edges) <- rownames(df1$edges)
  expect_that(df1$edges, equals(t1.e))

  expect_that(df1$vertices, equals(data.frame(row.names=letters[1:11],
                                              name=letters[1:11],
                                              stringsAsFactors=FALSE)))

  gg <- sg - g

  expect_that(ecount(gg), equals(0))
  expect_that(V(gg)$name, equals(letters[c(1:3,11)]))
})
