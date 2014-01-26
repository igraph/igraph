
#   IGraph R package
#   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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

# R function to plot compartment total curves from simul.net.epi .
# Inputs:  sim.res :=  list of simulated network SIR processes
#             comp := compartment (i.e., "NS", "NI", or "NR")
#                q := vector of lower and upper quantiles, resp
#             cols := char vector of colors for lines, median, and quantiles, resp.
# Outputs:  None.  Just produces the plot of all compartment curves, 
#           with median and quantiles.

plot.sir <- function(x, comp=c("NI", "NS", "NR"),
                     median=TRUE, quantiles=c(0.1, 0.9), color=NULL,
                     median_color=NULL, quantile_color=NULL,
                     lwd.median=2, lwd.quantile=2, lty.quantile=3,
                     xlim=NULL, ylim=NULL, xlab="Time", ylab=NULL, ...) {

  sir <- x

  if (!inherits(sir, "sir")) {
    stop("This is not an SIR model output")
  }
  comp <- toupper(igraph.match.arg(comp))
  if (!all(quantiles >= 0 & quantiles <= 1)) {
    stop("Quantiles should be in [0,1]")
  }
  
  if (is.null(color)) {
    color <- c(NI="skyblue", NS="pink", NR="palegoldenrod")[comp]
  }
  if (is.null(median_color)) {
    median_color <- c(NI="blue", NS="red", NR="gold")[comp]
  }
  if (is.null(quantile_color)) {
    quantile_color <- c(NI="blue", NS="red", NR="gold")[comp]
  }
  quantile_color <- rep(quantile_color, length.out=length(quantiles))

  ns <- length(sir)
  if (is.null(xlim)) {
    xlim <- c(0, max(sapply(sir, function(x) max(x$times))))
  }
  if (is.null(ylim)) {
    ylim <- c(0, max(sapply(sir, function(x) max(x[[comp]]))))
  }

  ## Work out median and quantile vectors.
  # Create single huge vector of times and compartment counts.
  if (median || length(quantiles) > 0) {
    big.time <- unlist(sapply(sir, function(x) { x$times }))
    big.N <- unlist(sapply(sir, function(x) { x[[comp]] }))
    ## Adhoc use of Freedman-Diaconis binwidth; rescale time accordingly.
    w <- 1/(2*(quantile(big.time, 0.75)-quantile(big.time, 0.25))/(100^(1/3)))
    time.bin <- floor(big.time*w)
  }
  
  ## Generate the plot, first with individual curves, and then 
  ## adding median and quantile curves.

  bin.vals <- as.numeric(levels(as.factor(time.bin)))/w

  if (is.null(ylab)) {
    if (comp == "NI") { ylab <- expression(N[I](t)) }
    if (comp == "NR") { ylab <- expression(N[R](t)) }
    if (comp == "NS") { ylab <- expression(N[S](t)) }
  }

  # Plot the stochastic curves individually.
  plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  lapply(seq_along(sir), function(i) {
    lines(sir[[i]]$time, sir[[i]][[comp]], col=color[1])
  })

  # Plot the median and quantiles.
  if (median) {
    my.m <- tapply(big.N, time.bin, median)
    lines(bin.vals, my.m, type="l", lwd=lwd.median, col=median_color)
  }
  for (i in seq_along(quantiles)) {
    my.ql <- tapply(big.N, time.bin, function(x) {
      quantile(x, prob=quantiles[i])
    })
    lines(bin.vals, my.ql, type="l", lty=lty.quantile, lwd=lwd.quantile,
          col=quantile_color[i])
  }

  invisible()
}
