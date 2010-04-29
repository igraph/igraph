
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
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

.onLoad <- function(dir, package) {
  library.dynam("igraph", package, dir, local=FALSE);

  ########################
  # Set default parameters
  ########################

  # printing attributes
  igraph.par("print.graph.attributes", FALSE)
  igraph.par("print.vertex.attributes", FALSE)
  igraph.par("print.edge.attributes", FALSE)

  # verbosity, progress bars mainly
  igraph.par("verbose", FALSE)

  # what to do with attributes
  igraph.par("vertex.attr.comb", list(name="concat", "ignore"))
  igraph.par("edge.attr.comb", list(weight="sum", name="concat", "ignore"))
}

.onUnload <- function(libpath) {
  library.dynam.unload("igraph", libpath)
}

.Last.lib <- function(libpath) {
  igraph::.onUnload(libpath)
}
