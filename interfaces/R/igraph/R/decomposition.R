#   IGraph R package
#   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

###################################################################
# Graph decomposition
###################################################################



#' Chordality of a graph
#' 
#' A graph is chordal (or triangulated) if each of its cycles of four or more
#' nodes has a chord, which is an edge joining two nodes that are not adjacent
#' in the cycle. An equivalent definition is that any chordless cycles have at
#' most three nodes.
#' 
#' The chordality of the graph is decided by first performing maximum
#' cardinality search on it (if the \code{alpha} and \code{alpham1} arguments
#' are \code{NULL}), and then calculating the set of fill-in edges.
#' 
#' The set of fill-in edges is empty if and only if the graph is chordal.
#' 
#' It is also true that adding the fill-in edges to the graph makes it chordal.
#'
#' @aliases is.chordal
#' @param graph The input graph. It may be directed, but edge directions are
#' ignored, as the algorithm is defined for undirected graphs.
#' @param alpha Numeric vector, the maximal chardinality ordering of the
#' vertices. If it is \code{NULL}, then it is automatically calculated by
#' calling \code{\link{max_cardinality}}, or from \code{alpham1} if
#' that is given..
#' @param alpham1 Numeric vector, the inverse of \code{alpha}. If it is
#' \code{NULL}, then it is automatically calculated by calling
#' \code{\link{max_cardinality}}, or from \code{alpha}.
#' @param fillin Logical scalar, whether to calculate the fill-in edges.
#' @param newgraph Logical scalar, whether to calculate the triangulated graph.
#' @return A list with three members: \item{chordal}{Logical scalar, it is
#' \code{TRUE} iff the input graph is chordal.} \item{fillin}{If requested,
#' then a numeric vector giving the fill-in edges. \code{NULL} otherwise.}
#' \item{newgraph}{If requested, then the triangulated graph, an \code{igraph}
#' object. \code{NULL} otherwise.}
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{max_cardinality}}
#' @references Robert E Tarjan and Mihalis Yannakakis. (1984). Simple
#' linear-time algorithms to test chordality of graphs, test acyclicity of
#' hypergraphs, and selectively reduce acyclic hypergraphs.  \emph{SIAM Journal
#' of Computation} 13, 566--579.
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## The examples from the Tarjan-Yannakakis paper
#' g1 <- graph_from_literal(A-B:C:I, B-A:C:D, C-A:B:E:H, D-B:E:F,
#'                 E-C:D:F:H, F-D:E:G, G-F:H, H-C:E:G:I,
#'                 I-A:H)
#' max_cardinality(g1)
#' is_chordal(g1, fillin=TRUE)
#' 
#' g2 <- graph_from_literal(A-B:E, B-A:E:F:D, C-E:D:G, D-B:F:E:C:G,
#'                 E-A:B:C:D:F, F-B:D:E, G-C:D:H:I, H-G:I:J,
#'                 I-G:H:J, J-H:I)
#' max_cardinality(g2)
#' is_chordal(g2, fillin=TRUE)
#' 
is_chordal <- function(graph, alpha = NULL, alpham1 = NULL,
                       fillin = FALSE, newgraph = FALSE) {
    if (!is_igraph(graph)) {
        stop("Not a graph object")
    }
    if (!is.null(alpha)) 
        alpha <- as.numeric(alpha)-1
    if (!is.null(alpham1)) 
        alpham1 <- as.numeric(alpham1)-1
    fillin <- as.logical(fillin)
    newgraph <- as.logical(newgraph)
    on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
    res <- .Call("R_igraph_is_chordal", graph, alpha, alpham1, 
                 fillin, newgraph, PACKAGE = "igraph")
    if (fillin) { res$fillin <- res$fillin + 1 }
    res
}
