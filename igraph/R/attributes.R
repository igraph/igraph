
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

add.graph.attribute <- function(graph, attrname) {
  .Call("R_igraph_add_graph_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

remove.graph.attribute <- function(graph, attrname) {
  .Call("R_igraph_remove_graph_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

get.graph.attribute <- function(graph, attrname=NULL) {
  if (is.null(attrname)) {
    res <- .Call("R_igraph_list_graph_attributes", graph,
                 PACKAGE="igraph")
    strsplit(rawToChar(res), "\n")[[1]]
  } else {
    .Call("R_igraph_get_graph_attribute", graph, as.character(attrname),
          PACKAGE="igraph")
  }
}

set.graph.attribute <- function(graph, attrname, value) {
  .Call("R_igraph_set_graph_attribute", graph, as.character(attrname),
        as.numeric(value),
        PACKAGE="igraph")
}

g.a <- get.graph.attribute
"g.a<-" <- set.graph.attribute

add.vertex.attribute <- function(graph, attrname) {
  .Call("R_igraph_add_vertex_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

remove.vertex.attribute <- function(graph, attrname) {
  .Call("R_igraph_remove_vertex_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

get.vertex.attribute <- function(graph, attrname=NULL,
                                 v=1:vcount(graph)-1) {
  if (is.null(attrname)) {
    res <- .Call("R_igraph_list_vertex_attributes", graph,
                 PACKAGE="igraph")
    strsplit(rawToChar(res), "\n")[[1]]
  } else {
    if (length(v)==1) {    
      .Call("R_igraph_get_vertex_attribute", graph,
            as.character(attrname), as.numeric(v),
            PACKAGE="igraph")
    } else {
      .Call("R_igraph_get_vertex_attributes", graph,
            as.character(attrname), as.numeric(v),
            PACKAGE="igraph")
    }
  }
}

set.vertex.attribute <- function(graph, attrname, v=1:vcount(graph)-1, value) {
  if (length(v)==1) {  
    .Call("R_igraph_set_vertex_attribute", graph, as.character(attrname),
          as.numeric(v), as.numeric(value),
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_set_vertex_attributes", graph, as.character(attrname),
          as.numeric(v), as.numeric(value),
          PACKAGE="igraph")
  }    
}

v.a <- get.vertex.attribute
"v.a<-" <- set.vertex.attribute

add.edge.attribute <- function(graph, attrname) {
  .Call("R_igraph_add_edge_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

remove.edge.attribute <- function(graph, attrname) {
  .Call("R_igraph_remove_edge_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

get.edge.attribute <- function(graph, attrname=NULL,
                                 e=1:ecount(graph)-1) {
  if (is.null(attrname)) {
    res <- .Call("R_igraph_list_edge_attributes", graph,
                 PACKAGE="igraph")
    strsplit(rawToChar(res), "\n")[[1]]
  } else {
    if (length(e)==1) {    
      .Call("R_igraph_get_edge_attribute", graph,
            as.character(attrname), as.numeric(e),
            PACKAGE="igraph")
    } else {
      .Call("R_igraph_get_edge_attributes", graph,
            as.character(attrname), as.numeric(e),
            PACKAGE="igraph")
    }
  }
}

set.edge.attribute <- function(graph, attrname, e=1:ecount(graph)-1, value) {
  if (length(e)==1) {  
    .Call("R_igraph_set_edge_attribute", graph, as.character(attrname),
          as.numeric(e), as.numeric(value),
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_set_edge_attributes", graph, as.character(attrname),
          as.numeric(e), as.numeric(value),
          PACKAGE="igraph")
  }    
}

e.a <- get.edge.attribute
"e.a<-" <- set.edge.attribute
