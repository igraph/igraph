
#   IGraph R package
#   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
#   Rue de l'Industrie 5, Lausanne 1005, Switzerland
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

.igraph.pb <- NULL

.igraph.progress <- function(percent, message) {
  type <- getIgraphOpt("verbose")
  if (is.logical(type) && type) {
    .igraph.progress.txt(percent, message)
  } else {
    switch (type,
            "tk"=.igraph.progress.tk(percent, message),
            stop("Cannot interpret 'verbose' option, this should not happen"))
  }
}

.igraph.progress.txt <- function(percent, message) {
  pb <- get(".igraph.pb", asNamespace("igraph"))
  if (percent==0) {
    if (!is.null(pb)) { close(pb) }
    cat(sep="", "  ", message, "\n")
    pb <- txtProgressBar(min=0, max=100, style=3)
  }
  setTxtProgressBar(pb, percent)
  if (percent==100) {
    close(pb);
    pb <- NULL
  }
  assign(".igraph.pb", pb, env=asNamespace("igraph"))
  0L
}

.igraph.progress.tk <- function(percent, message) {
  pb <- get(".igraph.pb", asNamespace("igraph"))
  if (percent==0) {
    if (!is.null(pb)) { close(pb) }
    pb <- tkProgressBar(min=0, max=100, title=message, label="0 %")
  }
  setTkProgressBar(pb, percent, label=paste(percent, "%"))
  if (percent==100) {
    close(pb);
    pb <- NULL
  }
  assign(".igraph.pb", pb, env=asNamespace("igraph"))
  0L
}
