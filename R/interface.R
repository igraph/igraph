
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

###################################################################
# Structure building
###################################################################

graph.empty <- function(..., type=igraph.par("default.type"))
  UseMethod(paste(sep="", "graph.empty.", type))

add.edges <- function(graph, edges)
  UseMethod(paste(sep="", "add.edges.", graph$gal$type))

add.vertices <- function(graph, nv)
  UseMethod(paste(sep="", "add.vertices.", graph$gal$type))

delete.edges <- function(graph, edges)
  UseMethod(paste(sep="", "delete.edges.", graph$gal$type))

delete.vertices <- function(graph, v)
  UseMethod(paste(sep="", "delete.vertices.", graph$gal$type))

###################################################################
# Structure query
###################################################################
  
vcount <- function(graph)
  UseMethod(paste(sep="", "vcount.", graph$gal$type))
  
ecount <- function(graph)
  UseMethod(paste(sep="", "ecount.", graph$gal$type))
  
neighbors <- function(graph, v, mode="out")
  UseMethod(paste(sep="", "neighbors.", graph$gal$type))

###################################################################
# Attributes, level 2
###################################################################
  
add.graph.attribute <- function(graph, attrname, default=NA)
  UseMethod(paste(sep="", "add.graph.attribute.", graph$gal$type))

delete.graph.attribute <- function(graph, attrname)
  UseMethod(paste(sep="", "delete.graph.attribute.", graph$gal$type))

get.graph.attribute <- function(graph, attrname=NULL)
  UseMethod(paste(sep="", "get.graph.attribute.", graph$gal$type))

set.graph.attribute <- function(graph, attrname, value)
  UseMethod(paste(sep="", "set.graph.attribute.", graph$gal$type))

add.vertex.attribute <- function(graph, attrname, type="simple", default=NA)
  UseMethod(paste(sep="", "add.vertex.attribute.", graph$gal$type))

delete.vertex.attribute <- function(graph, attrname)
  UseMethod(paste(sep="", "delete.vertex.attribute.", graph$gal$type))

get.vertex.attribute <- function(graph, attrname=NULL, v=NULL)
  UseMethod(paste(sep="", "get.vertex.attribute.", graph$gal$type))

set.vertex.attribute <- function(graph, attrname, v=NULL, value)
  UseMethod(paste(sep="", "set.vertex.attribute.", graph$gal$type))

add.edge.attribute <- function(graph, attrname, type="simple", default=NA)
  UseMethod(paste(sep="", "add.edge.attribute.", graph$gal$type))

delete.edge.attribute <- function(graph, attrname)
  UseMethod(paste(sep="", "delete.edge.attribute.", graph$gal$type))

get.edge.attribute <- function(graph, attrname=NULL, from=NULL, to=NULL)
  UseMethod(paste(sep="", "get.edge.attribute.", graph$gal$type))

set.edge.attribute <- function(graph, attrname, from=NULL, to=NULL, value)
  UseMethod(paste(sep="", "set.edge.attribute.", graph$gal$type))

g.a <- get.graph.attribute
"g.a<-" <- set.graph.attribute
v.a <- get.vertex.attribute
"v.a<-" <- set.vertex.attribute
e.a <- get.edge.attribute
"e.a<-" <- set.edge.attribute

###################################################################
# Iterators, level 2
###################################################################

igraph.iterator <- function(graph, type="vid")
  UseMethod(paste(sep="", "igraph.iterator.", graph$gal$type))

###################################################################
# Interface object for C functions
###################################################################

igraph.c.interface <- do.call
