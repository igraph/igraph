## -----------------------------------------------------------------------
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
## -----------------------------------------------------------------------


#' Palette for categories
#'
#' This is a color blind friendly palette from
#' \url{http://jfly.iam.u-tokyo.ac.jp/color}. It has 8 colors.
#'
#' This is the suggested palette for visualizations where vertex colors
#' mark categories, e.g. community membership.
#'
#' @param n The number of colors in the palette. We simply take the first
#' \code{n} colors from the total 8.
#' @return A character vector of RGB color codes.
#'
#' @family palettes
#' @export
#' @examples
#' library(igraphdata)
#' data(karate)
#' karate <- karate %>%
#'   add_layout_(with_fr()) %>%
#'   set_vertex_attr("size", value = 10)
#'
#' cl_k <- cluster_optimal(karate)
#'
#' V(karate)$color <- membership(cl_k)
#' karate$palette <- categorical_pal(length(cl_k))
#' plot(karate)

categorical_pal <- function(n) {

  stopifnot(n > 0)

  x <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
         "#D55E00", "#CC79A7", "#999999")

  if (n > length(x)) warning("Cannot make ", n, " categorical colors")

  n <- min(n, length(x))

  x[seq_len(n)]
}


#' Sequential palette
#'
#' This is the \sQuote{OrRd} palette from \url{http://colorbrewer2.org}.
#' It has at most nine colors.
#'
#' Use this palette, if vertex colors mark some ordinal quantity, e.g. some
#' centrality measure, or some ordinal vertex covariate, like the age of
#' people, or their seniority level.
#'
#' @param n The number of colors in the palette. The maximum is nine
#' currently.
#' @return A character vector of RGB color codes.
#'
#' @family palettes
#' @export
#' @examples
#' library(igraphdata)
#' data(karate)
#' karate <- karate %>%
#'   add_layout_(with_kk()) %>%
#'   set_vertex_attr("size", value = 10)
#'
#' V(karate)$color <- scales::dscale(degree(karate) %>% cut(5), sequential_pal)
#' plot(karate)

sequential_pal <- function(n) {

  stopifnot(n >= 0)

  x <- list(
    "#FEE8C8",
    c("#FEE8C8", "#FDBB84"),
    c("#FEE8C8", "#FDBB84", "#E34A33"),
    c("#FEF0D9", "#FDCC8A", "#FC8D59", "#D7301F"),
    c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"),
    c("#FEF0D9", "#FDD49E", "#FDBB84", "#FC8D59", "#E34A33", "#B30000"),
    c("#FEF0D9", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F",
      "#990000"),
    c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548",
      "#D7301F", "#990000"),
    c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548",
      "#D7301F", "#B30000", "#7F0000")
  )

  if (n > length(x)) warning("Cannot make ", n, " sequential colors")

  n <- min(n, length(x))

  if (n == 0) character() else x[[n]]
}


#' Diverging palette
#'
#' This is the \sQuote{PuOr} palette from \url{http://colorbrewer2.org}.
#' It has at most eleven colors.
#'
#' This is similar to \code{\link{sequential_pal}}, but it also puts
#' emphasis on the mid-range values, plus the the two extreme ends.
#' Use this palette, if you have such a quantity to mark with vertex
#' colors.
#'
#' @param n The number of colors in the palette. The maximum is eleven
#' currently.
#' @return A character vector of RGB color codes.
#'
#' @family palettes
#' @export
#' @examples
#' library(igraphdata)
#' data(foodwebs)
#' fw <- foodwebs[[1]] %>%
#'   induced_subgraph(V(.)[ECO == 1]) %>%
#'   add_layout_(with_fr()) %>%
#'   set_vertex_attr("label", value = seq_len(gorder(.))) %>%
#'   set_vertex_attr("size", value = 10) %>%
#'   set_edge_attr("arrow.size", value = 0.3)
#'
#' V(fw)$color <- scales::dscale(V(fw)$Biomass %>% cut(10), diverging_pal)
#' plot(fw)
#'
#' library(igraphdata)
#' data(karate)
#' karate <- karate %>%
#'   add_layout_(with_kk()) %>%
#'   set_vertex_attr("size", value = 10)
#'
#' V(karate)$color <- scales::dscale(degree(karate) %>% cut(5), diverging_pal)
#' plot(karate)

diverging_pal <- function(n) {

  stopifnot(n > 0)

  x <- list(
    "#F1A340",
    c("#F1A340", "#F7F7F7"),
    c("#F1A340", "#F7F7F7", "#998EC3"),
    c("#E66101", "#FDB863", "#B2ABD2", "#5E3C99"),
    c("#E66101", "#FDB863", "#F7F7F7", "#B2ABD2", "#5E3C99"),
    c("#B35806", "#F1A340", "#FEE0B6", "#D8DAEB", "#998EC3", "#542788"),
    c("#B35806", "#F1A340", "#FEE0B6", "#F7F7F7", "#D8DAEB", "#998EC3",
      "#542788"),
    c("#B35806", "#E08214", "#FDB863", "#FEE0B6", "#D8DAEB", "#B2ABD2",
      "#8073AC", "#542788"),
    c("#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7", "#D8DAEB",
      "#B2ABD2", "#8073AC", "#542788"),
    c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#D8DAEB",
      "#B2ABD2", "#8073AC", "#542788", "#2D004B"),
    c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7",
      "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")
  )

  if (n > length(x)) warning("Cannot make ", n, " divergent colors")

  n <- min(n, length(x))

  if (n == 0) character() else x[[n]]
}


#' The default R palette
#'
#' This is the default R palette, to be able to reproduce the
#' colors of older igraph versions. Its colors are appropriate
#' for categories, but they are not very attractive.
#'
#' @param n The number of colors to use, the maximum is eight.
#' @return A character vector of color names.
#'
#' @family palettes
#' @export

r_pal <- function(n) {
  x <- palette()
  if (n > length(x)) warning("Cannot make ", n, " divergent colors")
  n <- min(n, length(x))
  if (n == 0) character() else x[[n]]
}
