\name{plot.sir}
\alias{plot.sir}
\concept{SIR model}
\title{Plotting the results on multiple SIR model runs}
\description{
  This function can conveniently plot the results of multiple SIR model
  simulations.
}
\usage{
\S3method{plot}{sir}(x, comp = c("NI", "NS", "NR"), median = TRUE,
    quantiles = c(0.1, 0.9), color = NULL, median_color = NULL,
    quantile_color = NULL, lwd.median = 2, lwd.quantile = 2,
    lty.quantile = 3, xlim = NULL, ylim = NULL,
    xlab = "Time", ylab = NULL, \dots)
}
\arguments{
  \item{x}{The output of the SIR simulation, coming from the
    \code{\link{sir}} function.}
  \item{comp}{Character scalar, which component to plot. Either
    \sQuote{NI} (infected, default), \sQuote{NS} (susceptible) or
    \sQuote{NR} (recovered). }
  \item{median}{Logical scalar, whether to plot the (binned) median.}
  \item{quantiles}{A vector of (binned) quantiles to plot.}
  \item{color}{Color of the individual simulation curves.}
  \item{median_color}{Color of the median curve.}
  \item{quantile_color}{Color(s) of the quantile curves. (It is recycled
    if needed and non-needed entries are ignored if too long.)}
  \item{lwd.median}{Line width of the median.}
  \item{lwd.quantile}{Line width of the quantile curves.}
  \item{lty.quantile}{Line type of the quantile curves.}
  \item{xlim}{The x limits, a two-element numeric vector. If
    \code{NULL}, then it is calculated from the data.}
  \item{ylim}{The y limits, a two-element numeric vector. If
    \code{NULL}, then it is calculated from the data.}
  \item{xlab}{The x label.}
  \item{ylab}{The y label. If \code{NULL} then it is automatically added
    based on the \code{comp} argument.}
  \item{\dots}{Additional arguments are passed to \code{plot}, that is
    run before any of the curves are added, to create the figure.}
}
\details{
  The number of susceptible/infected/recovered individuals is
  plotted over time, for multiple simulations.
}
\value{Nothing.}
\references{
  Bailey, Norman T. J. (1975). The mathematical theory of infectious
  diseases and its applications (2nd ed.). London: Griffin.
}
\author{
  Eric Kolaczyk (\url{http://math.bu.edu/people/kolaczyk/}) and
  Gabor Csardi \email{csardi.gabor@gmail.com}.
}
\seealso{\code{\link{sir}} for running the actual simulation.}
\examples{
g <- erdos.renyi.game(100, 100, type="gnm")
sm <- sir(g, beta=5, gamma=1)
plot(sm)
}
\keyword{graphs}
  