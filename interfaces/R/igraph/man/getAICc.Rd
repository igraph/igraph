\name{getAICc}
\alias{getAICc}
\concept{Compute Information Criteria Value Using a Poisson Approximation}
\title{Compute AIC based on a Poisson Approximation using the output from \code{gclust}}
\description{
    Compute and Extract information Criteria Value from \code{gclust} using a Poisson
    approximation, where the penality term is adjusted for small sample cases.  
}
\usage{ getAICc(gfit) }
\arguments{
    \item{gfit}{The output from a call to \code{gclust.rsvt} or \code{gclust.app}}
}

\value{
    \item{nclust}{Number of clusters being considered}
    \item{negloglikpart}{Negative log likelihood part}
    \item{parampart}{Parameter penalty term}
    \item{AIC}{AIC Value Using a Possion approximation}
}

\author{ Nam Lee \email{nhlee@jhu.edu}}

\examples{
require(NMF)
require(irlba)
require(igraph)
}

