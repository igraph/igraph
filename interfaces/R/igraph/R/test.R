#   IGraph R package
#   Copyright (C) 2005-2013  Gabor Csardi <csardi.gabor@gmail.com>
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



#' Run package tests
#' 
#' Runs all package tests.
#' 
#' The \code{testthat} package is needed to run all tests. The location tests
#' themselves can be extracted from the package via \code{system.file("tests",
#' package="igraph")}.
#' 
#' This function simply calls the \code{test_dir} function from the
#' \code{testthat} package on the test directory.
#'
#' @aliases igraphtest
#' @return Whatever is returned by \code{test_dir} from the \code{testthat}
#' package.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
#' @export

igraph_test <- function() {
  do.call(require, list("testthat"))
  tdir <- system.file("tests", package="igraph")
  do.call("test_dir", list(tdir))
}



#' Query igraph's version string
#' 
#' Queries igraph's original version string. See details below.
#' 
#' The igraph version string is the same as the version of the R package for
#' all realeased igraph versions. For development versions and nightly builds,
#' they might differ however.
#' 
#' The reason for this is, that R package version numbers are not flexible
#' enough to cover in-between releases versions, e.g. alpha and beta versions,
#' release candidates, etc.
#'
#' @aliases igraph.version
#' @return A character scalar, the igraph version string.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
#' @export
#' @examples
#' 
#' ## Compare to the package version
#' packageDescription("igraph")$Version
#' igraph_version()

igraph_version <- function() {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_version", PACKAGE="igraph")
}

checkpkg <- function(package_file, args=character()) {
  package_file <- as.character(package_file)
  args <- as.character(args)
  do.call(":::", list("tools", ".check_packages"))(c(package_file, args))
}
