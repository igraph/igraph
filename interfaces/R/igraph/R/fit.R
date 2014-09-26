#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

###################################################################
# Pit a power-law (khmm a Yule really) distribution,
# this is a common degree distribution in networks
###################################################################



#' Fitting a power-law distribution function to discrete data
#' 
#' \code{power.law.fit} fits a power-law distribution to a data set.
#' 
#' This function fits a power-law distribution to a vector containing samples
#' from a distribution (that is assumed to follow a power-law of course). In a
#' power-law distribution, it is generally assumed that \eqn{P(X=x)} is
#' proportional to \eqn{x^{-alpha}}{x^-alpha}, where \eqn{x} is a positive
#' number and \eqn{\alpha}{alpha} is greater than 1. In many real-world cases,
#' the power-law behaviour kicks in only above a threshold value
#' \eqn{x_{min}}{xmin}. The goal of this function is to determine
#' \eqn{\alpha}{alpha} if \eqn{x_{min}}{xmin} is given, or to determine
#' \eqn{x_{min}}{xmin} and the corresponding value of \eqn{\alpha}{alpha}.
#' 
#' \code{power.law.fit} provides two maximum likelihood implementations.  If
#' the \code{implementation} argument is \sQuote{\code{R.mle}}, then the BFGS
#' optimization (see \link[stats4]{mle}) algorithm is applied.  The additional
#' arguments are passed to the mle function, so it is possible to change the
#' optimization method and/or its parameters.  This implementation can
#' \emph{not} to fit the \eqn{x_{min}}{xmin} argument, so use the
#' \sQuote{\code{plfit}} implementation if you want to do that.
#' 
#' The \sQuote{\code{plfit}} implementation also uses the maximum likelihood
#' principle to determine \eqn{\alpha}{alpha} for a given \eqn{x_{min}}{xmin};
#' When \eqn{x_{min}}{xmin} is not given in advance, the algorithm will attempt
#' to find itsoptimal value for which the \eqn{p}-value of a Kolmogorov-Smirnov
#' test between the fitted distribution and the original sample is the largest.
#' The function uses the method of Clauset, Shalizi and Newman to calculate the
#' parameters of the fitted distribution. See references below for the details.
#' 
#' @param x The data to fit, a numeric vector. For implementation
#' \sQuote{\code{R.mle}} the data must be integer values. For the
#' \sQuote{\code{plfit}} implementation non-integer values might be present and
#' then a continuous power-law distribution is fitted.
#' @param xmin Numeric scalar, or \code{NULL}. The lower bound for fitting the
#' power-law. If \code{NULL}, the smallest value in \code{x} will be used for
#' the \sQuote{\code{R.mle}} implementation, and its value will be
#' automatically determined for the \sQuote{\code{plfit}} implementation. This
#' argument makes it possible to fit only the tail of the distribution.
#' @param start Numeric scalar. The initial value of the exponent for the
#' minimizing function, for the \sQuote{\code{R.mle}} implementation. Ususally
#' it is safe to leave this untouched.
#' @param force.continuous Logical scalar. Whether to force a continuous
#' distribution for the \sQuote{\code{plfit}} implementation, even if the
#' sample vector contains integer values only (by chance). If this argument is
#' false, igraph will assume a continuous distribution if at least one sample
#' is non-integer and assume a discrete distribution otherwise.
#' @param implementation Character scalar. Which implementation to use. See
#' details below.
#' @param \dots Additional arguments, passed to the maximum likelihood
#' optimizing function, \code{\link[stats4]{mle}}, if the \sQuote{\code{R.mle}}
#' implementation is chosen. It is ignored by the \sQuote{\code{plfit}}
#' implementation.
#' @return Depends on the \code{implementation} argument. If it is
#' \sQuote{\code{R.mle}}, then an object with class \sQuote{\code{mle}}. It can
#' be used to calculate confidence intervals and log-likelihood. See
#' \code{\link[stats4]{mle-class}} for details.
#' 
#' If \code{implementation} is \sQuote{\code{plfit}}, then the result is a
#' named list with entries: \item{continuous}{Logical scalar, whether the
#' fitted power-law distribution was continuous or discrete.}
#' \item{alpha}{Numeric scalar, the exponent of the fitted power-law
#' distribution.} \item{xmin}{Numeric scalar, the minimum value from which the
#' power-law distribution was fitted. In other words, only the values larger
#' than \code{xmin} were used from the input vector.} \item{logLik}{Numeric
#' scalar, the log-likelihood of the fitted parameters.} \item{KS.stat}{Numeric
#' scalar, the test statistic of a Kolmogorov-Smirnov test that compares the
#' fitted distribution with the input vector. Smaller scores denote better
#' fit.} \item{KS.p}{Numeric scalar, the p-value of the Kolmogorov-Smirnov
#' test. Small p-values (less than 0.05) indicate that the test rejected the
#' hypothesis that the original data could have been drawn from the fitted
#' power-law distribution.}
#' @author Tamas Nepusz \email{ntamas@@gmail.com} and Gabor Csardi
#' \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link[stats4]{mle}}
#' @references Power laws, Pareto distributions and Zipf's law, M. E. J.
#' Newman, \emph{Contemporary Physics}, 46, 323-351, 2005.
#' 
#' Aaron Clauset, Cosma R .Shalizi and Mark E.J. Newman: Power-law
#' distributions in empirical data. SIAM Review 51(4):661-703, 2009.
#' @keywords graphs
#' @examples
#' 
#' # This should approximately yield the correct exponent 3
#' g <- barabasi.game(1000) # increase this number to have a better estimate
#' d <- degree(g, mode="in")
#' fit1 <- power.law.fit(d+1, 10)
#' fit2 <- power.law.fit(d+1, 10, implementation="R.mle")
#' 
#' fit1$alpha
#' coef(fit2)
#' fit1$logLik
#' logLik(fit2)
#' 
power.law.fit <- function(x, xmin=NULL, start=2, force.continuous=FALSE,
                          implementation=c("plfit", "R.mle"), ...) {

  implementation <- igraph.match.arg(implementation)

  if (implementation == "r.mle") {
    power.law.fit.old(x, xmin, start, ...)
  } else if (implementation == "plfit") {
    if (is.null(xmin)) xmin <- -1
    power.law.fit.new(x, xmin=xmin, force.continuous=force.continuous)
  }
}

power.law.fit.old <- function(x, xmin=NULL, start=2, ...) {

  if (length(x) == 0) {
    stop("zero length vector")
  }
  if (length(x) == 1) {
    stop("vector should be at least of length two")
  }  

  require(stats4)
  
  if (is.null(xmin)) { xmin <- min(x) }
  
  n <- length(x)
  x <- x[ x >= xmin]
  if (length(x) != n) {
    n <- length(x)
  }
  
#  mlogl <- function(alpha) {
#    if (xmin > 1) {
#      C <- 1/(1/(alpha-1)-sum(beta(1:(xmin-1), alpha)))
#    } else {
#      C <- alpha-1
#    }
#    -n*log(C)-sum(lbeta(x, alpha))
#  }

  mlogl <- function(alpha) {
     C <- 1/sum( (xmin:10000)^-alpha )
     -n*log(C)+alpha*sum(log(x))
  }

  alpha <- mle(mlogl, start=list(alpha=start), ...)

  alpha
}
