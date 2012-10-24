
###############################################
## Part 1, initialization, don't run it yet

library(websockets)
library(igraph)

w <- createContext(port=7681L, webpage=NULL)

f <- function(DATA, WS, ...) {
  if (is.raw(DATA)) { DATA <- rawToChar(DATA); }
  print(paste("Received:", DATA))
  send_graph(mygraph, WS)
}

set_callback('receive', f, w)

send_graph <- function(graph, ws) {
  if (!is.igraph(graph)) {
    stop("Not an igraph graph")
  }

  ## Vertices
  if (is.named(graph)) {
    ids <- V(graph)$name
  } else {
    ids <- as.character(seq_len(vcount(graph)))
  }
  msg <- sprintf('{"an":{"%s":{"label":"%s"}}}', ids, ids)
  fmsg1 <- paste(msg, collapse="\n\r")

  ## Edges
  el <- get.edgelist(graph, names=TRUE)
  mode(el) <- "character"
  msg <- sprintf('{"ae":{"%s:%s":{"source":"%s","target":"%s"}}}',
                 el[,1], el[,2], el[,1], el[,2])
  fmsg2 <- paste(msg, collapse="\n\r")
  fmsg <- paste(fmsg1, sep="\n\r", fmsg2)
  websocket_write(fmsg, ws)
}

cl <- function(WS) {
  print("Closed")
}
set_callback('closed', cl, w)

es <- function(WS) {
  print("Established")
}

set_callback('established',es,w)

add.node <- function() {
  msg <- sprintf('{"an":{"%d":{"label":"%d"}}}',
                 vcount(mygraph)+1, vcount(mygraph)+1)
  mygraph <<- mygraph + 1
  websocket_broadcast(msg, w)
}

add.edge <- function(from, to) {
  msg <- sprintf('{"ae":{"%d:%d":{"source":"%d","target":"%d"}}}',
                 from, to, from, to)
  mygraph <<- mygraph + edge(from, to)
  websocket_broadcast(msg, w)
}

########################################
## Part 2, set the initial graph, and start running the server

# mygraph <- nexus.get("karate")
mygraph <- graph(c(1,2,2,3))

daemonize(w)

########################################
## Part 3, make some operations, from R

add.node()
add.edge(1,4)

########################################
## Part 4, stop

websocket_close(w)
