
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

add.graph.attribute <- function(graph, attrname, type="numeric") {
  if (is.character(type)) {
    type <- switch(type, "numeric"=0, "character"=1, "string"=1)
  }
  .Call("R_igraph_add_graph_attribute", graph, as.character(attrname),
        as.numeric(type), PACKAGE="igraph")
}

remove.graph.attribute <- function(graph, attrname) {
  .Call("R_igraph_remove_graph_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

get.graph.attribute <- function(graph, attrname=NULL) {
  if (is.null(attrname)) {
    .Call("R_igraph_list_graph_attributes", graph,
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_get_graph_attribute", graph, as.character(attrname),
          PACKAGE="igraph")
  }
}

set.graph.attribute <- function(graph, attrname, value) {
  .Call("R_igraph_set_graph_attribute", graph, as.character(attrname),
        value, PACKAGE="igraph")
}

g.a <- get.graph.attribute
"g.a<-" <- set.graph.attribute

add.vertex.attribute <- function(graph, attrname, type="numeric") {
  if (is.character(type)) {
    type <- switch(type, "numeric"=0, "character"=1, "string"=1)
  }
  .Call("R_igraph_add_vertex_attribute", graph, as.character(attrname),
        as.numeric(type), PACKAGE="igraph")
}

remove.vertex.attribute <- function(graph, attrname) {
  .Call("R_igraph_remove_vertex_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

get.vertex.attribute <- function(graph, attrname=NULL,
                                 v=igraph.vs.all(graph)) {
  if (is.null(attrname)) {
    .Call("R_igraph_list_vertex_attributes", graph,
          PACKAGE="igraph")
  } else {
    if (is.numeric(v) && length(v)==1) {    
      .Call("R_igraph_get_vertex_attribute", graph,
            as.character(attrname), as.numeric(v),
            PACKAGE="igraph")
    } else {
      .Call("R_igraph_get_vertex_attributes", graph,
            as.character(attrname), as.igraph.vs(graph, v),
            PACKAGE="igraph")
    }
  }
}

set.vertex.attribute <- function(graph, attrname,
                                 v=igraph.vs.all(graph), value) {
  if (is.numeric(v) && length(v)==1) {  
    .Call("R_igraph_set_vertex_attribute", graph, as.character(attrname),
          as.numeric(v), value,
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_set_vertex_attributes", graph, as.character(attrname),
          as.igraph.vs(graph, v), value,
          PACKAGE="igraph")
  }    
}

v.a <- get.vertex.attribute
"v.a<-" <- set.vertex.attribute

add.edge.attribute <- function(graph, attrname, type="numeric") {
  if (is.character(type)) {
    type <- switch(type, "numeric"=0, "character"=1, "string"=1)
  }
  .Call("R_igraph_add_edge_attribute", graph, as.character(attrname),
        as.numeric(type), PACKAGE="igraph")
}

remove.edge.attribute <- function(graph, attrname) {
  .Call("R_igraph_remove_edge_attribute", graph, as.character(attrname),
        PACKAGE="igraph")
}

get.edge.attribute <- function(graph, attrname=NULL,
                                 e=igraph.es.all(graph)) {
  if (is.null(attrname)) {
    .Call("R_igraph_list_edge_attributes", graph,
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_get_edge_attributes", graph,
          as.character(attrname), e,
          PACKAGE="igraph")
  }
}

set.edge.attribute <- function(graph, attrname, e=igraph.es.all(graph),value) {
  .Call("R_igraph_set_edge_attributes", graph, as.character(attrname),
        e, value,
        PACKAGE="igraph")
}

e.a <- get.edge.attribute
"e.a<-" <- set.edge.attribute
