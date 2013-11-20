
benchmark_context(
  id="local.scan",
  version="1"
  )

benchmark_this("",
  init = {
    library(igraph)
    set.seed(42)
  },
  init_each = {
    g <- random.graph.game(1000, p=.2)    
  },
  replications = 10,
  scan0 = {
    ls0 <- local.scan(g, k=0)
    error
  },
  scan1 = {
    ls1 <- local.scan(g)
  })
