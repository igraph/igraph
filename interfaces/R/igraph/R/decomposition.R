
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

is.chordal <- function(graph, alpha = NULL, alpham1 = NULL,
                       fillin = FALSE, newgraph = FALSE) {
    if (!is.igraph(graph)) {
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
