
#   IGraph R package
#   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

###################################################################
# Common functions for plot and tkplot
###################################################################

i.get.layout <- function(graph, layout, layout.par) {

  if (is.function(layout)) {
    layout <- layout(graph, layout.par)
  } else if (is.character(layout) && length(layout)==1 &&
             substr(layout, 1, 2)=="a:") {
    layout <- matrix(unlist(get.vertex.attribute(graph, substring(layout,3))),
                     nr=vcount(graph), byrow=TRUE)[,1:2]
  }  

  layout
}

i.get.vertex.color <- function(graph, vertex.color) {

  if (length(vertex.color)==1 && substr(vertex.color, 1, 2)=="a:") {
    vertex.color <- unlist(get.vertex.attribute(graph,
                                                substring(vertex.color,3)))
  }

  if (is.numeric(vertex.color)) {
    vertex.color <- vertex.color %% length(palette())
    vertex.color[vertex.color==0] <- length(palette())
    vertex.color <- palette()[vertex.color]
  }

  vertex.color  
}

i.get.vertex.frame.color <- function(graph, vertex.frame.color) {

  if (length(vertex.frame.color)==1 &&
      substr(vertex.frame.color, 1, 2)=="a:") {
    vertex.frame.color <-
      unlist(get.vertex.attribute(graph,
                                  substring(vertex.frame.color,3)))
  }

  if (is.numeric(vertex.frame.color)) {
    vertex.frame.color <- vertex.frame.color %% length(palette())
    vertex.frame.color[vertex.frame.color==0] <- length(palette())
    vertex.frame.color <- palette()[vertex.frame.color]
  }

  vertex.frame.color  
}

i.get.vertex.size <- function(graph, vertex.size) {

  if (is.character(vertex.size) &&
      length(vertex.size)==1 && substr(vertex.size, 1, 2)=="a:") {
    vertex.size <- as.numeric(get.vertex.attribute
                                 (graph, substring(vertex.size,3)))
  }
  vertex.size
}

i.get.edge.color <- function(graph, edge.color) {

  if (length(edge.color)==1 && substr(edge.color, 1, 2)=="a:") {
    edge.color <- as.character(get.edge.attribute
                               (graph, substring(edge.color,3)))
  }

  if (is.numeric(edge.color)) {
    edge.color <- edge.color %% length(palette())
    edge.color[edge.color==0] <- length(palette())
    edge.color <- palette()[edge.color]
  }
  edge.color
}

i.get.edge.width <- function(graph, edge.width) {

  if (is.character(edge.width) &&
      length(edge.width)==1 && substr(edge.width, 1, 2)=="a:") {
    edge.width <- as.character(get.edge.attribute
                               (graph, substring(edge.width,3)))
  }
  edge.width
}

i.get.edge.labels <- function(graph, edge.labels) {

  if (is.character(edge.labels) &&
      length(edge.labels)==1 && substr(edge.labels, 1, 2)=="a:") {
    edge.labels <- as.character(get.edge.attribute
                               (graph, substring(edge.labels,3)))
  }
  edge.labels
}

i.get.label.degree <- function(graph, label.degree) {

  if (is.character(label.degree) &&
      length(label.degree)==1 && substr(label.degree, 1, 2)=="a:") {
    label.degree <- as.numeric(get.vertex.attribute
                               (graph, substring(label.degree,3)))
  }
  label.degree
}

i.get.labels <- function(graph, labels) {

  if (is.null(labels)) {
    labels <- 0:(vcount(graph)-1)
  } else if (is.na(labels[1])) {
  } else if (is.character(labels) && length(labels)==1 &&
             substr(labels, 1, 2)=="a:") {
    labels <- unlist(get.vertex.attribute(graph, substring(labels, 3)))
  }
  labels
}
