
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

ba.game <- function(n, m=NULL, out.dist=NULL, out.seq=NULL,
                    out.pref=FALSE, ...) {

  # Checks
  if (! is.null(out.seq) && (!is.null(m) || !is.null(out.dist))) {
    warning("if `out.seq' is given `m' and `out.dist' should be NULL")
    m <- out.dist <- NULL
  }
  if (is.null(out.seq) && !is.null(out.dist) && !is.null(m)) {
    warning("if `out.dist' is given `m' will be ignored")
    m <- NULL
  }
  if (!is.null(out.seq) && length(out.seq) != n) {
    stop("`out.seq' should be of length `n'")
  }
  if (!is.null(out.seq) && min(out.seq)<0) {
    stop("negative elements in `out.seq'");
  }
  if (!is.null(m) && m<0) {
    stop("`m' is negative")
  }
  if (!is.null(m) && m==0) {
    warning("`m' is zero, graph will be empty")
  }

  if (is.null(m) && is.null(out.dist) && is.null(out.seq)) {
    m <- 1
  }
  
  n <- as.numeric(n)
  if (!is.null(m)) { m <- as.numeric(m) }
  if (!is.null(out.dist)) { out.dist <- as.numeric(out.dist) }
  if (!is.null(out.dist)) { out.seq <- as.numeric(out.seq) }
  out.pref <- as.logical(out.pref)

  if (!is.null(out.dist)) {
    out.seq <- as.numeric(sample(0:(length(out.dist)-1), n,
                                 replace=TRUE, prob=out.dist))
  }
  
  res <- graph.empty(n=n, ...)
  res <- add.edges(res, .Call("REST_ba_game", n, m, out.seq, out.pref,
                              PACKAGE="igraph"))
  
  res
}

barabasi.game <- ba.game
