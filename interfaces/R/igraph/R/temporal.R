
#   IGraph R package
#   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
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
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

as.time <- function(graph, time, calendar=graph$calendar, NA.OK=FALSE) {
  res <- ifelse(time == Inf, -1L, match(time, calendar))
  if (!NA.OK && any(is.na(res))) {
    stop("Unknown time point")
  }
  res
}

check_calendar <- function(calendar) {
  if (any(duplicated(calendar))) {
    warning("Making calendar unique")
    calendar <- unique(calendar)
  }
  if (is.unsorted(calendar)) {
    warning("Sorting calendar")
    calendar <- sort(calendar)
  }
  calendar
}

is.temporal <- function(graph) {
  stopifnot(is.igraph(graph))
  "calendar" %in% list.graph.attributes(graph)
}

as.temporal <- function(graph, calendar = 0L) {
  calendar <- check_calendar(calendar)
  stopifnot(is.igraph(graph))
  graph$calendar <- calendar
  graph
}

graph.temporal <- function(edges, n = max(edges), directed = TRUE,
                          e_on = 0, e_off = NULL, v_on = 0, v_off = NULL,
                          calendar = unique(sort(c(e_on, e_off, v_on,
                            v_off)))) {

  ec <- length(edges) / 2
  stopifnot(length(e_on) == ec || length(e_on) == 1)
  stopifnot(length(v_on) == n  || length(v_on) == 1)
  stopifnot(is.null(e_off) || length(e_off) == ec || length(e_off) == 1)
  stopifnot(is.null(v_off) || length(v_off) == n  || length(v_off) == 1)

  if (!missing(n) && n < max(edges)) { n <- max(edges) }

  if (!missing(calendar)) { calendar <- check_calendar(calendar) }

  check <- function(v, len) {
    if (!is.null(v)) {
      v <- match(v, calendar)
      if (any(is.na(v))) { stop("Unknown time label") }
      v <- rep(v, length.out=len)
    }
    v
  }

  e_on <- check(e_on, ec)
  e_off <- check(e_off, ec)
  v_on <- check(v_on, n)
  v_off <- check(v_off, n)

  on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
  res <- .Call("R_igraph_create_temporal", as.numeric(edges) - 1,
              as.numeric(n), as.logical(directed), e_on, e_off, v_on,
              v_off, PACKAGE="igraph")
  res$calendar <- calendar

  res
}

on <- function(set) {
  isvs <- inherits(set, "igraph.vs")
  ises <- inherits(set, "igraph.es")
  if (!isvs && !ises) {
    stop("Not a set of vertices or a set of edges")
  }

  graph <- attr(set, "env")$graph

  res <- if (isvs) {
    graph.vertices.range(graph, set, on=TRUE, off=FALSE)
  } else {
    graph.edges.range(graph, set, on=TRUE, off=FALSE)
  }

  res[ res == -1L ] <- NA
  ifelse(is.na(res), Inf, graph$calendar[res])
}

off <- function(set) {
  isvs <- inherits(set, "igraph.vs")
  ises <- inherits(set, "igraph.es")
  if (!isvs && !ises) {
    stop("Not a set of vertices or a set of edges")
  }

  graph <- attr(set, "env")$graph

  res <- if (isvs) {
    graph.vertices.range(graph, set, on=FALSE, off=TRUE)
  } else {
    graph.edges.range(graph, set, on=FALSE, off=TRUE)
  }

  res[ res == -1L ] <- NA
  ifelse(is.na(res), Inf, graph$calendar[res])
}

graph.vertices.range <- function(graph, vids=V(graph), on=TRUE, off=TRUE) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)
  if (!on && !off) {
    stop("At least on of 'on' and 'off' must be TRUE")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_vertices_range", graph, vids-1, as.logical(on),
              as.logical(off), PACKAGE="igraph")

  res
}

graph.edges.range <- function(graph, eids=E(graph), on=TRUE, off=TRUE) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  eids <- as.igraph.es(graph, eids)
  if (!on && !off) {
    stop("At least on of 'on' and 'off' must be TRUE")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_edges_range", graph, eids-1, as.logical(on),
              as.logical(off), PACKAGE="igraph")

  res
}
