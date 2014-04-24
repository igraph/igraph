
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

igraphtest <- function() {
  do.call(require, list("testthat"))
  tdir <- system.file("tests", package="igraph")
  do.call("test_dir", list(tdir))
}

igraph.version <- function() {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_version", PACKAGE="igraph")
}

checkpkg <- function(package_file, args=character()) {
  package_file <- as.character(package_file)
  args <- as.character(args)
  do.call(":::", list("tools", ".check_packages"))(c(package_file, args))
}
