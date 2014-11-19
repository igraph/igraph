
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

make_call <- function(f, ..., .args = list()) {
  if (is.character(f)) f <- as.name(f)
  as.call(c(f, ..., .args))
}

do_call <- function(f, ..., .args = list(), .env = parent.frame()) {
  f <- substitute(f)

  call <- make_call(f, ..., .args)
  eval(call, .env)
}

add_class <- function(x, class) {
  if (!is(x, class)) {
    class(x) <- c(class, class(x))
  }
  x
}

`%||%` <- function (lhs, rhs) {
  lres <- withVisible(eval(lhs, envir = parent.frame()))
  if (is.null(lres$value)) {
    eval(rhs, envir = parent.frame())
  } else {
    if (lres$visible) {
      lres$value
    } else {
      invisible(lres$value)
    }
  }
}

`%&&%` <- function(lhs, rhs) {
  lres <- withVisible(eval(lhs, envir = parent.frame()))
  if (!is.null(lres$value)) {
    eval(rhs, envir = parent.frame())
  } else {
    if (lres$visible) {
      lres$value
    } else {
      invisible(lres$value)
    }
  }
}

## Grab all arguments of the parent call, in a list

grab_args <- function() {
  envir <- parent.frame()
  func <- sys.function(-1)
  call <- sys.call(-1)
  dots <- match.call(func, call, expand.dots=FALSE)$...
  c(as.list(envir), dots)
}

capitalize <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

address <- function(x) {
  .Call("R_igraph_address", x, PACKAGE = "igraph")
}

`%+%` <- function(x, y) {
  stopifnot(is.character(x), is.character(y))
  paste0(x, y)
}

chr <- as.character

drop_null <- function(x) {
  x [!sapply(x, is.null)]
}
