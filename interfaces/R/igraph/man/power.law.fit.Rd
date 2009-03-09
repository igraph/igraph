\name{power.law.fit}
\alias{power.law.fit}
\concept{Power-law}
\title{Fitting a power-law distribution function to discrete data}
\description{\code{power.law.fit} fits a power-law distribution to a
  data set.}
\usage{
power.law.fit(x, xmin = NULL, start = 2, \dots)
}
\arguments{
  \item{x}{The data to fit, a numeric vector containing integer values.}
  \item{xmin}{The lower bound for fitting the power-law. If NULL, the
    smallest value in \code{x} will be used. This argument makes it
    possible to fit only the tail of the distribution. }
  \item{start}{The initial value of the exponent for the minimizing
    function. Ususally it is safe to leave this untouched.}
  \item{\dots}{Additional arguments, passed to the maximum likelyhood
    optimizing function, \code{\link[stats4]{mle}}.}
}
\details{
  A power-law distribution is fitted with maximum likelyhood
  methods as recommended by Newman and (by default) the 
  BFGS optimization (see \link[stats4]{mle}) algorithm is applied.

  The additional arguments are passed to the mle function, so it is
  possible to change the optimization method and/or its parameters.
}
\value{An object with class \sQuote{\code{mle}}. It can be used to
  calculate confidence intervals and log-likelihood. See
  \code{\link[stats4]{mle-class}} for details.}
\references{Power laws, Pareto distributions and Zipf's law,
  M. E. J. Newman, \emph{Contemporary Physics}, in press.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link[stats4]{mle}}}
\examples{
# This should approximately yield the correct exponent 3
g <- barabasi.game(1000)    # increase this number to have a better estimation
d <- degree(g, mode="in")
power.law.fit(d+1, 20)
}
\keyword{graphs}
