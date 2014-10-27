
#' @useDynLib igraph
#' @import methods
#' @importFrom magrittr %>%
#' @export make_bipartite_graph
#' @export connect
#' @export make_de_bruijn_graph
#' @export make_full_bipartite_graph
#' @export graph_from_adjacency_matrix
#' @export graph_from_data_frame
#' @export graph_from_incidence_matrix
#' @export is_matching
#' @export is_max_matching
#' @export make_kautz_graph
#' @export make_line_graph
#' @export max_bipartite_match
#' @export sample_asym_pref
NULL

#' Magrittr's pipes
#'
#' igraph re-exports the \code{\%>\%} operator of magrittr, because
#' we find it very useful. Please see the documentation in the
#' \code{magrittr} package.
#'
#' @param lhs Left hand side of the pipe.
#' @param rhs Right hand side of the pipe.
#' @return Result of applying the right hand side to the
#'   result of the left hand side.
#'
#' @export
#' @rdname pipe
#' @examples
#' make_ring(10) %>%
#'   add_edges(c(1,6)) %>%
#'   plot()

`%>%` <- magrittr::`%>%`
