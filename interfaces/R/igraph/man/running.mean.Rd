\name{running.mean}
\alias{running.mean}
\title{Running mean of a time series}
\description{\code{running.mean} calculates the running mean in a vector
  with the given bin width. The running mean is defined by the mean
  of the previous \code{binwidth} values.}
\usage{
running.mean(v, binwidth)
}
\arguments{
  \item{v}{The numeric vector.}
  \item{binwidth}{Numeric constant, the size of the bin.}
}
% \details{}
\value{
  A numeric vector of length \code{length(v)-binwidth+1}
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
running.mean(1:1000, 10)
}
\keyword{manip}
