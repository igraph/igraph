
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

# TODO: clean disconnect

igraph.par("pgsql.host", "localhost")
igraph.par("pgsql.dbname", "igraph")
igraph.par("pgsql.user", "igraph")
igraph.par("pgsql.user", NULL)

###################################################################
# Structure building
###################################################################

graph.empty.pgsql.default <- function(..., type="pgsql",
                                      host=igraph.par("pgsql.host"),
                                      dbname=igraph.par("pgsql.dbname"),
                                      user=igraph.par("pgsql.user")) {

  require(Rdbi.PgSQL)

  gal <- list(type=type, ...)
  if (! "directed" %in% names(gal)) {
    gal$directed <- TRUE
  }
  realn <- ifelse(is.null(gal$n), 0, gal$n)
  gal$n <- 0

  # connect only if not connected
  conn <- connect.if.needed(host, dbname, user)
  
  dbSendQuery(conn, "BEGIN")
  qres <- dbGetQuery(conn, "select nextval('graph_gname_seq')")
  gname <- as.numeric(qres[1])
  attributes(gname) <- NULL
  dbSendQuery(conn, "INSERT INTO graph (gname) VALUES (", gname, ")")
  
  for (i in seq(along=gal)) {
    dbSendQuery(conn, "INSERT INTO gal (gname, name, value) VALUES (",
                gname, ",", sQuote(names(gal)[i]), ",",
                sQuote(serialize(gal[[i]], connection=NULL, ascii=TRUE)), ")")
  }
  dbSendQuery(conn, "COMMIT")

  data <- list(gname=gname, conn=conn, host=host, dbname=dbname, user=user)
  res <- list(data=data, gal=list(type=type))
  class(res) <- "graph"
  res <- add.vertices(res, realn)

  res
}

add.edges.pgsql.default <- function(graph, edges) {

  add.edges.common(graph, edges)

  edges <- data.frame(graph$data$gname, matrix(edges, nc=2, byrow=TRUE))
  names(edges) <- c("gname", "v1", "v2")
  dbAppendTable(graph$data$conn, "data", edges)
  
  graph
}

add.vertices.pgsql.default <- function(graph, nv) {

  add.vertices.common(graph, nv)

  n <- vcount(graph)
  set.graph.attribute(graph, "n", n+nv)
  
  graph
}

delete.edges.pgsql.default <- function(graph, edges) {

  delete.edges.common(graph, edges)

  dbSendQuery(graph$data$conn, "BEGIN")
  # TODO: do this with one query only (possible?)
  idx <- 1
  if (is.directed(graph)) {
    while (idx <= length(edges)) {
      dbSendQuery(graph$data$conn, "DELETE from data where gname=",
                  graph$data$gname, 
                  "and eid=(select min(eid) from data where gname=",
                  graph$data$gname, "and v1=", edges[idx],
                  "and v2=", edges[idx+1], ")" )
      idx <- idx + 2
    }
  } else {
    while (idx <= length(edges)) {
      dbSendQuery(graph$data$conn, "DELETE from data where gname=",
                  graph$data$gname, 
                  "and eid=(select min(eid) from data where gname=",
                  graph$data$gname, "and ((v1=", edges[idx],
                  "and v2=", edges[idx+1], ") or (v1=", edges[idx+1],
                  "and v2=", edges[idx], ")))")
      idx <- idx + 2
    }
  }
  dbSendQuery(graph$data$conn, "COMMIT")

  graph
}

delete.vertices.pgsql.default <- function(graph, v) {

  v <- unique(v)
  nodes <- vcount(g)
  replment <- rev(rev((1:nodes)[-v])[1:length(v)])
  print (v)
  print (replment)

  dbSendQuery(graph$data$conn, "BEGIN")
  for (i in v) {
    dbSendQuery(graph$data$conn, "DELETE from data where gname=",
                graph$data$gname, "and v1=", i, "or v2=", i)
  }
  for (i in seq(along=v)) {
    dbSendQuery(graph$data$conn, "UPDATE data set v1=", v[i],
                "where gname=", graph$data$gname, "and v1=", replment[i]) 
    dbSendQuery(graph$data$conn, "UPDATE data set v2=", v[i],
                "where gname=", graph$data$gname, "and v2=", replment[i])
  }
  dbSendQuery(graph$data$conn, "COMMIT")

  graph
}

###################################################################
# Structure query
###################################################################

vcount.pgsql.default <- function(graph) {
  get.graph.attribute(graph, "n")
}
    
ecount.pgsql.default <- function(graph) {

  qres <- dbGetQuery(graph$data$conn, "select count(*) from data where gname=",
                     graph$data$gname)
  res <- qres[[1]]
  attributes(res) <- NULL

  res
}
    
neighbors.pgsql.default <- function(graph, v, mode="out") {

  res <- numeric()
  dbSendQuery(graph$data$conn, "BEGIN")  
  if (!is.directed(graph) || mode=="all" || mode=="out") {
    qres <- dbGetQuery(graph$data$conn, "select v2 from data where gname=",
                       graph$data$gname, "and v1=", v)
    res <- c(res, qres$v2)
  }
  if (!is.directed(graph) || mode=="all" || mode=="in") {
    qres <- dbGetQuery(graph$data$conn, "select v1 from data where gname=",
                       graph$data$gname, "and v2=", v)
    res <- c(res, qres$v1)
  }
  dbSendQuery(graph$data$conn, "COMMIT")

  res
}

###################################################################
# Attributes
###################################################################
  
add.graph.attribute.pgsql.default <- function(graph, attrname, default=NA) {

  dbSendQuery(graph$data$conn,
              "insert into gal (gname, name, value) values (",
              graph$data$gname, ",", sQuote(attrname), ",",
              sQuote(serialize(default, connection=NULL, ascii=TRUE)),
              ")")
  
  graph
}

delete.graph.attribute.pgsql.default <- function(graph, attrname) {

  dbSendQuery(graph$data$conn, "delete from gal where gname=",
              graph$data$gname, "and name=", sQuote(attrname))

  graph
}

get.graph.attribute.pgsql.default <- function(graph, attrname=NULL) {

  if (is.null(attrname)) {
    qres <- dbGetQuery(graph$data$conn,
                       "select distinct name from gal where gname=",
                       graph$data$gname)
    res <- qres[,"name"]
  } else {
    qres <- dbGetQuery(graph$data$conn, "select value from gal where gname=",
                       graph$data$gname, "and name=", sQuote(attrname))
    if (length(qres$value)==0) {
      res <- NULL
    } else {
      res <- unserialize(as.character(qres[1]))
    }
  }
  
  res
}

set.graph.attribute.pgsql.default <- function(graph, attrname, value) {

  dbSendQuery(graph$data$conn, "update gal set value=",
              sQuote(serialize(value, connection=NULL, ascii=TRUE)),
              "where gname=", graph$data$gname,
              "and name=", sQuote(attrname))

  graph
}

add.vertex.attribute.pgsql.default <- function(graph, attrname, default=NA) {

  n <- vcount(graph)
  value <- sQuote(serialize(default, connection=NULL, ascii=TRUE))
  dbSendQuery(graph$data$conn, "BEGIN")
  dbSendQuery(graph$data$conn,
              "insert into val (gname, v, name, value) values (",
              graph$data$gname, ",", 0, ",", sQuote(attrname), ",",
              value, ")")
  
  if (n>=1) {
    for (i in 1:n) {
      dbSendQuery(graph$data$conn,
                  "insert into val (gname, v, name, value) values (",
                  graph$data$gname, ",", i, ",", sQuote(attrname), ",",
                  value, ")")
    }
  }
  dbSendQuery(graph$data$conn, "COMMIT")
  
  graph
}


delete.vertex.attribute.pgsql.default <- function(graph, attrname) {

  dbSendQuery(graph$data$conn, "delete from val where gname=",
              graph$data$gname, "and name=", sQuote(attrname))

  graph
}

get.vertex.attribute.pgsql.default <- function(graph, attrname=NULL, v=NULL) {

  if (is.null(attrname)) {
    qres <- dbGetQuery(graph$data$conn,
                       "select distinct name from val where gname=",
                       graph$data$gname)
    res <- qres[,"name"]
  } else {
    dbSendQuery(graph$data$conn, "BEGIN")
    if (is.null(v)) {
      qres <- dbGetQuery(graph$data$conn,
                         "select value from val where gname=",
                         graph$data$gname,
                         "and name=", sQuote(attrname),
                         "and v<>0 order by v")
      res <- unname(sapply(qres$value, unserialize, simplify=FALSE))
    } else {
      res <- vector(length(v), mode="list")
      for (i in seq(along=v)) {
        qres <- dbGetQuery(graph$data$conn,
                           "select value from val where gname=",
                           graph$data$gname, "and name=", sQuote(attrname),
                           "and v=", v[i])
        res[[i]] <- unname(unserialize(qres$value))
        rm(qres)
      }
    }
    dbSendQuery(graph$data$conn, "COMMIT")
    
    if (length(res)==0) {
      stop("no such attribute: ", attrname)
    }
    
    if (length(res)==1) {
      res <- res[[1]]
    }
  }

  res
}

set.vertex.attribute.pgsql.default <- function(graph, attrname, v=NULL,
                                               value) {
  
  dbSendQuery(graph$data$conn, "BEGIN")
  if (is.null(v)) {
    dbSendQuery(graph$data$conn, "update val set value=",
                sQuote(serialize(value, connection=NULL, ascii=TRUE)),
                "where gname=", graph$data$gname, "and name=",
                sQuote(attrname))
  } else if (length(v) != 0) {
    dbSendQuery(graph$data$conn, "update val set value=",
                sQuote(serialize(value, connection=NULL, ascii=TRUE)),
                "where gname=", graph$data$gname, "and name=",
                sQuote(attrname), "and v in (",
                paste(v, collapse=","), ")")
  }
  dbSendQuery(graph$data$conn, "COMMIT")
  
  graph  
}

# TODO
add.edge.attribute.pgsql.default <- function(graph, attrname, default=NA) {
  error("This storage type does not support edge attributes (yet)")
}

# TODO
delete.edge.attribute.pgsql.default <- function(graph, attrname) {
  error("This storage type does not support edge attributes (yet)")
}

# TODO
get.edge.attribute.pgsql.default <- function(graph, attrname=NULL,
                                             from=NULL, to=NULL) {
  error("This storage type does not support edge attributes (yet)")
}

# TODO
set.edge.attribute.pgsql.default <- function(graph, attrname,
                                             from=NULL, to=NULL, value) {
  error("This storage type does not support edge attributes (yet)")
}

###################################################################
# Storage type specific functions
###################################################################

# Lists all graphs in a database
pgsql.list.graphs <- function(host=igraph.par("pgsql.host"),
                              dbname=igraph.par("pgsql.dbname"),
                              user=igraph.par("pgsql.user")) {

  require(Rdbi.PgSQL)

  conn <- connect.if.needed(host, dbname, user)
  gnames <- dbGetQuery(conn, "SELECT gname FROM graph")
  unname(gnames$gname)
}

# Returns a graph with the given gname, ... for future use
pgsql.get.graph <- function(gname,
                            host=igraph.par("pgsql.host"),
                            dbname=igraph.par("pgsql.dbname"),
                            user=igraph.par("pgsql.user"), ...) {
  
  require(Rdbi.PgSQL)

  conn <- connect.if.needed(host, dbname, user)
  gnames <- dbGetQuery(conn, "SELECT gname FROM graph")
  if (!gname %in% gnames$gname) {
    warning("Could not find graph in database, proceeding")
  }
  res <- list(data=list(conn=conn, gname=gname, host=host,
                 dbname=dbname, user=user), gal=list(type="pgsql"))
  class(res) <- "graph"
  
  res
}

# Deletes all temporary graphs, or all graphs from a database
pgsql.clear <- function(all=FALSE,
                        host=igraph.par("pgsql.host"),
                        dbname=igraph.par("pgsql.dbname"),
                        user=igraph.par("pgsql.user")) {
  
  require(Rdbi.PgSQL)

  conn <- connect.if.needed(host, dbname, user)

  dbSendQuery(conn, "BEGIN")
  if (all) {
    dbSendQuery(conn, "delete from graph")
  } else {
    graphs <- pgsql.list.graphs(host=host, dbname=dbname, user=user)
    for (gname in graphs) {
      g <- pgsql.get.graph(gname=gname, host=host, dbname=dbname, user=user)
      temp <- g.a(g, "temporary")
      if (!is.null(temp) && is.logical(temp) && temp) {
        dbSendQuery("delete from graph where gname=", gname)
      }
    }
  }
  dbSendQuery(conn, "COMMIT")

  invisible()
}

###################################################################
# Internal functions
###################################################################

connect.if.needed <- function(host, dbname, user) {
  
  connect <- TRUE
  conn <- igraph.par("pgsql.conn")
  if (!is.null(conn) && inherits(conn, "PgSQL.conn")) {
    conninfo <- dbConnectionInfo(conn)
    if (conninfo$database.name == dbname &&
        conninfo$host.name == host &&
        conninfo$user.name == user &&
        conninfo$status == 0) {
      connect <- FALSE
    } 
  }
  if (connect) {
    conn <- dbConnect.PgSQL(host=host, dbname=dbname, user=user)
    igraph.par("pgsql.conn", conn)
  }
  
  return(conn);
}
