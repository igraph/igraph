
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
  
  g <- graph.disjoint.union(g1, g2)

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

test_that("disjoint union gives warning for non-unique vertex names", {

  library(igraph)

  g1 <- graph.ring(5); V(g1)$name <- letters[1:5]
  g2 <- graph.ring(5); V(g2)$name <- letters[5:9]
  
  expect_that(graph.disjoint.union(g1, g2),
              gives_warning("Duplicate vertex names in disjoint union"))
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

  g <- graph.intersection(g1, g2, keep.all.vertices=FALSE)

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

test_that("compose works for named graphs", {

  library(igraph)

  g1 <- graph.formula( A-B:D:E, B-C:D, C-D, D-E )
  g2 <- graph.formula( A-B-E-A )

  V(g1)$bar1 <- seq_len(vcount(g1))
  V(g2)$bar2 <- seq_len(vcount(g2))
  V(g1)$foo <- letters[seq_len(vcount(g1))]
  V(g2)$foo <- letters[seq_len(vcount(g2))]

  E(g1)$bar1 <- seq_len(ecount(g1))
  E(g2)$bar2 <- seq_len(ecount(g2))
  E(g1)$foo <- letters[seq_len(ecount(g1))]
  E(g2)$foo <- letters[seq_len(ecount(g2))]

  g <- graph.compose(g1, g2)
  df <- get.data.frame(g, what="both")

  df.v <- read.table(stringsAsFactors=FALSE, textConnection("
  bar1 foo_1 foo_2 bar2 name
A    1     a     a    1    A
B    2     b     b    2    B
D    3     c    NA   NA    D
E    4     d     c    3    E
C    5     e    NA   NA    C
"))
  expect_that(df$vertices, equals(df.v))

  df.e <- read.table(stringsAsFactors=FALSE, textConnection("
   from to bar1 foo_1 foo_2 bar2
1     A  B    3     c     c    3
2     A  A    3     c     b    2
3     A  E    1     a     c    3
4     A  A    1     a     a    1
5     B  E    1     a     b    2
6     B  B    1     a     a    1
7     B  D    6     f     c    3
8     A  D    6     f     b    2
9     D  E    4     d     c    3
10    A  D    4     d     a    1
11    D  E    2     b     b    2
12    B  D    2     b     a    1
13    E  E    3     c     b    2
14    B  E    3     c     a    1
15    E  C    5     e     c    3
16    A  C    5     e     a    1
"))
  rownames(df$edges) <- rownames(df$edges)
  expect_that(df$edges, equals(df.e))

})

test_that("intersection of non-named graphs keeps attributes properly", {
  library(igraph)
  set.seed(42)

  g <- erdos.renyi.game(10, 1/2)
  g2 <- erdos.renyi.game(10, 1/2)
  E(g)$weight <- sample(ecount(g))
  E(g2)$weight <- sample(ecount(g2))

  gi <- graph.intersection(g, g2)

  rn <- function(D) {
    rownames(D) <- paste(D[,1], D[,2], sep="-")
    D
  }

  df <- rn(get.data.frame(g))
  df2 <- rn(get.data.frame(g2))
  dfi <- rn(get.data.frame(gi))

  expect_that(df[rownames(dfi), ], is_equivalent_to(dfi[, 1:3]))
  expect_that(df2[rownames(dfi), ], is_equivalent_to(dfi[, c(1,2,4)]))

})

test_that("union of non-named graphs keeps attributes properly", {
  library(igraph)
  set.seed(42)

  g <- erdos.renyi.game(10, 1/2)
  g2 <- erdos.renyi.game(10, 1/2)
  E(g)$weight <- sample(ecount(g))
  E(g2)$weight <- sample(ecount(g2))

  gu <- graph.union(g, g2)

  rn <- function(D) {
    rownames(D) <- paste(D[,1], D[,2], sep="-")
    D
  }

  df <- rn(get.data.frame(g))
  df2 <- rn(get.data.frame(g2))
  dfu <- rn(get.data.frame(gu))

  expect_that(dfu[rownames(df), 1:3], is_equivalent_to(df))
  expect_that(dfu[rownames(df2), c(1,2,4)], is_equivalent_to(df2))

  expect_that(dfu[!rownames(dfu) %in% rownames(df), 3],
              equals(rep(NA_real_, ecount(gu)-ecount(g))))
  expect_that(dfu[!rownames(dfu) %in% rownames(df2), 4],
              equals(rep(NA_real_, ecount(gu)-ecount(g2))))

})
