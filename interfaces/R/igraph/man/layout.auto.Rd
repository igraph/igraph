\name{layout.auto}
\alias{layout.auto}
\concept{Graph layout}
\title{Choose an appropriate graph layout algorithm automaticall}
\description{This function tries to choose an appropriate
  graph layout algorithm for the graph, automatically, based on a simple
  algorithm. See details below.}
\usage{
layout.auto (graph, dim = 2, \dots)
}
\arguments{
  \item{graph}{The input graph}
  \item{dim}{Dimensions, should be 2 or 3.}
  \item{\dots}{The extra arguments are passed to the real layout
    function.}
}
\details{
  \code{layout.auto} tries to choose an appropriate layout function for
  the supplied graph, and uses that to generate the layout. The current
  implementation works like this:
  \enumerate{
    \item If the graph has a graph attribute called \sQuote{layout},
    then this is used. If this attribute is an R function, then it is
    called, with the graph and any other extra arguments.
    \item Otherwise, if the graph has vertex attributes called
    \sQuote{x} and \sQuote{y}, then these are used as coordinates. If
    the graph has an additional \sQuote{z} vertex attribute, that is
    also used.
    \item Otherwise, if the graph is connected and has less than 1000
    vertices, the Kamada-Kawai layout is used, by calling
    \code{layout.kamada.kawai}.
    \item Otherwise the DrL layout is used, \code{layout.drl} is called.
  } 
}
\value{
  A numeric matrix with two or three columns.
}
% \references{}
\seealso{\code{\link{plot.igraph}}}
% \examples{}
\keyword{graphs}
