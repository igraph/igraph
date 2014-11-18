
## ----------------------------------------------------------------------
##
##   IGraph R package
##   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
##   334 Harvard street, Cambridge, MA 02139 USA
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##   02110-1301 USA
##
## ----------------------------------------------------------------------

#' Igraph data structure versions
#'
#' Igraph's internal data representation changes sometimes between
#' versions. This means that it is not possible to use igraph objects
#' that were created (and possibly saved to a file) with an older
#' igraph version.
#'
#' \code{graph_version} queries the current data format,
#' or the data format of a possibly older igraph graph.
#'
#' \code{\link{upgrade_graph}} can convert an older data format
#' to the current one.
#'
#' @param graph The input graph. If it is missing, then
#'   the version number of the current data format is returned.
#' @return A character scalar.
#'
#' @seealso upgrade_graph to convert the data format of a graph.
#' @export

graph_version <- function(graph) {
  if (missing(graph)) {
    "0.8.0"

  } else {
    stopifnot(is_igraph(graph))
    .Call("R_igraph_graph_version", graph, PACKAGE = "igraph")
  }
}

#' Igraph data structure versions
#'
#' Igraph's internal data representation changes sometimes between
#' versions. This means that it is not possible to use igraph objects
#' that were created (and possibly saved to a file) with an older
#' igraph version.
#'
#' \code{\link{graph_version}} queries the current data format,
#' or the data format of a possibly older igraph graph.
#'
#' \code{upgrade_graph} can convert an older data format
#' to the current one.
#'
#' @param graph The input graph.
#' @return The graph in the current format.
#'
#' @seealso graph_version to check the current data format version
#' or the version of a graph.
#' @export

upgrade_graph <- function(graph) {

  stopifnot(is_igraph(graph))

  g_ver <- graph_version(graph)
  p_ver <- graph_version()

  if (g_ver < p_ver) {

    if (g_ver == "0.4.0" && p_ver == "0.8.0") {
      .Call("R_igraph_add_env", graph, PACKAGE = "igraph")

    } else {
      stop("Don't know how to upgrade graph from ", g_ver, " to ", p_ver)
    }

  } else if (g_ver > p_ver) {
    stop("Don't know how to downgrade graph from ", g_ver, " to ", p_ver)

  } else {
    graph
  }
}
