#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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



#' Running mean of a time series
#' 
#' \code{running.mean} calculates the running mean in a vector with the given
#' bin width.
#' 
#' The running mean of \code{v} is a \code{w} vector of length
#' \code{length(v)-binwidth+1}. The first element of \code{w} id the average of
#' the first \code{binwidth} elements of \code{v}, the second element of
#' \code{w} is the average of elements \code{2:(binwidth+1)}, etc.
#' 
#' @param v The numeric vector.
#' @param binwidth Numeric constant, the size of the bin, should be meaningful,
#' ie. smaller than the length of \code{v}.
#' @return A numeric vector of length \code{length(v)-binwidth+1}
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords manip
#' @examples
#' 
#' running.mean(1:100, 10)
#' 
running.mean <- function(v, binwidth) {

  v <- as.numeric(v)
  binwidth <- as.numeric(binwidth)
  if (length(v) < binwidth) {
    stop("Vector too short for this binwidth.")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_running_mean", v, binwidth,
       PACKAGE="igraph");
}



#' Sampling a random integer sequence
#' 
#' This function provides a very efficient way to pull an integer random sample
#' sequence from an integer interval.
#' 
#' The algorithm runs in \code{O(length)} expected time, even if
#' \code{high-low} is big. It is much faster (but of course less general) than
#' the builtin \code{sample} function of R.
#' 
#' @param low The lower limit of the interval (inclusive).
#' @param high The higher limit of the interval (inclusive).
#' @param length The length of the sample.
#' @return An increasing numeric vector containing integers, the sample.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @references Jeffrey Scott Vitter: An Efficient Algorithm for Sequential
#' Random Sampling, \emph{ACM Transactions on Mathematical Software}, 13/1,
#' 58--67.
#' @keywords datagen
#' @examples
#' 
#' rs <- igraph.sample(1, 100000000, 10)
#' rs
#' 
igraph.sample <- function(low, high, length) {
  if (length>high-low+1) {
    stop("length too big for this interval")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_random_sample", as.numeric(low), as.numeric(high),
        as.numeric(length),
        PACKAGE="igraph")
}

igraph.match.arg <- function(arg, choices, several.ok=FALSE) {
  if (missing(choices)) {
    formal.args <- formals(sys.function(sys.parent()))
    choices <- eval(formal.args[[deparse(substitute(arg))]])
  }

  arg <- tolower(arg)
  choices <- tolower(choices)

  match.arg(arg=arg, choices=choices, several.ok=several.ok)
}

igraph.i.spMatrix <- function(M) {
  if (M$type == "triplet") {
    Matrix::sparseMatrix(dims=M$dim, i=M$i+1L, j=M$p+1L, x=M$x)
  } else {
    new("dgCMatrix", Dim=M$dim, Dimnames=list(NULL, NULL),
        factors=list(), i=M$i, p=M$p, x=M$x)
  }
}



#' Set random seed of the C library's RNG
#' 
#' Set the random seed of the C library's RNG, for a new sequence of
#' pseudo-random numbers.
#' 
#' Note that this function has nothing to do with R's random number generator,
#' see \code{set.seed} for that.
#' 
#' Some package (e.g. ngspatial) use internal C code and generate random
#' numbers using the standard C library's built-in random number generator
#' instead of using R's RNGs. The \code{srand} function is provided to set the
#' random seed for these packages. It simply calls the standard C function
#' \code{srand}, with the supplied integer seed value.
#' 
#' Note that the standard C library's RNGs are typically of very bad quality,
#' and also slower than R's RNGs. It is not worth using them, really, other
#' than taking over some legacy C code that already uses them, and that would
#' be difficult to rewrite to use R's RNGs.
#' 
#' @param seed Numeric scalar, the new random seed. It must be non-negative and
#' will be converted to an integer.
#' @return \code{NULL}, invisibly.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
srand <- function(seed) {
  seed <- as.numeric(seed)
  if (length(seed) != 1) { stop("Length of `seed' must be 1") }
  if (seed < 0) { stop("Seed must be non-negative") }
  res <- .Call("R_igraph_srand", seed, PACKAGE="igraph")
  invisible(res)
}
