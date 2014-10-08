
#   IGraph R package
#   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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

# This is a sparse data frame. It is like a regular data frame,
# but it allows for some columns to be constant, and then it
# stores that column more economically.

sdf <- function(..., row.names = NULL, NROW = NULL) {
  
  cols <- list(...)
  
  if (is.null(names(cols)) || any(names(cols) == "") ||
      any(duplicated(names(cols)))) {
    stop("Columns must be have (unique) names")
  }
  
  lens <- sapply(cols, length)
  n1lens <- lens[ lens != 1 ]
  
  if (length(unique(n1lens)) > 1) {
    stop("Columns must be constants or have the same length")
  }
  
  if (length(n1lens) == 0) {
    if (is.null(NROW)) {
      stop("Cannot determine number of rows")
    }
    attr(cols, "NROW") <- NROW
  } else {
    if (!is.null(NROW) && n1lens[1] != NROW) {
      stop("NROW does not match column lengths")
    }
    attr(cols, "NROW") <- unname(n1lens[1])
  }

  class(cols) <- "igraphSDF"
  attr(cols, "row.names") <- row.names
  
  cols
}

#' @method as.data.frame igraphSDF

as.data.frame.igraphSDF <- function(x, row.names, optional, ...) {
  as.data.frame(lapply(x, rep, length.out=attr(x, "NROW")))
}

#' @method "[" igraphSDF

`[.igraphSDF` <- function(x, i, j, ..., drop=TRUE) {
  if (!is.character(j)) {
    stop("The column index must be character")
  }
  if (!missing(i) && !is.numeric(i)) {
    stop("The row index must be numeric")
  }
  if (missing(i)) {
    rep(x[[j]], length.out=attr(x, "NROW"))
  } else {
    if (length(x[[j]])==1) {
      rep(x[[j]], length(i))
    } else {
      x[[j]][i]
    }
  }
}

#' @method "[<-" igraphSDF
 
`[<-.igraphSDF` <- function(x, i, j, value) {
  if (!is.character(j)) {
    stop("The column index must be character")
  }
  if (!missing(i) && !is.numeric(i)) {
    stop("Row index must be numeric, if given")
  }
  if (missing(i)) {
    if (length(value) != attr(x, "NROW") && length(value) != 1) {
      stop("Replacement value has the wrong length")
    }
    x[[j]] <- value
  } else {
    if (length(value) != length(i) && length(value) != 1) {
      stop("Replacement value has the wrong length")
    }
    tmp <- rep(x[[j]], length=attr(x, "NROW"))
    tmp[i] <- value
    if (length(unique(tmp)) == 1) {
      tmp <- tmp[1]
    }
    x[[j]] <- tmp
  }

  x
}
