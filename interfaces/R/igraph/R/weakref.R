
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

make_weak_ref <- function(key, value, finalizer = NULL) {
  base::.Call("R_igraph_make_weak_ref", key, value,
              finalizer, PACKAGE = "igraph")
}

weak_ref_key <- function(ref) {
  base::.Call("R_igraph_weak_ref_key", ref, PACKAGE = "igraph")
}

weak_ref_value <- function(ref) {
  base::.Call("R_igraph_weak_ref_value", ref, PACKAGE = "igraph")
}

weak_ref_run_finalizer <- function(ref) {
  base::.Call("R_igraph_weak_ref_run_finalizer", ref,
              PACKAGE = "igraph")
}
