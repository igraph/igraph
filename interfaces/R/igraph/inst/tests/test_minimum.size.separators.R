
context("minimum.size.separators")

test_that("minimum.size.separators works", {

  library(igraph)

  camp <- graph.formula(Harry:Steve:Don:Bert - Harry:Steve:Don:Bert,
                        Pam:Brazey:Carol:Pat - Pam:Brazey:Carol:Pat,
                        Holly   - Carol:Pat:Pam:Jennie:Bill,
                        Bill    - Pauline:Michael:Lee:Holly,
                        Pauline - Bill:Jennie:Ann,
                        Jennie  - Holly:Michael:Lee:Ann:Pauline,
                        Michael - Bill:Jennie:Ann:Lee:John,
                        Ann     - Michael:Jennie:Pauline,
                        Lee     - Michael:Bill:Jennie,
                        Gery    - Pat:Steve:Russ:John,
                        Russ    - Steve:Bert:Gery:John,
                        John    - Gery:Russ:Michael)
  camp <- simplify(camp)
  sep <- lapply(minimum.size.separators(camp), function(x) V(camp)[x])
  expect_that(all(sapply(sep, is.minimal.separator, graph=camp)), is_true())

})
