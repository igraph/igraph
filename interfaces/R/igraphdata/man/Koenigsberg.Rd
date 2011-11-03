\name{Koenigsberg}
\docType{data}
\alias{Koenigsberg}
\title{Bridges of Koenigsberg from Euler's times}
\description{
  The Seven Bridges of Koenigsberg is a notable historical problem in
  mathematics. Its negative resolution by Leonhard Euler in 1735 laid
  the foundations of graph theory and presaged the idea of topology. 

  The city of Koenigsberg in Prussia (now Kaliningrad, Russia) was set on
  both sides of the Pregel River, and included two large islands which
  were connected to each other and the mainland by seven bridges

  The problem was to find a walk through the city that would cross each
  bridge once and only once. The islands could not be reached by any route
  other than the bridges, and every bridge must have been crossed
  completely every time (one could not walk half way onto the bridge and
  then turn around and later cross the other half from the other side).

  Euler proved that the problem has no solution.
}
\usage{
Koenigsberg
}
\format{
  An undirected \code{igraph} graph object with vertex attributes
  \sQuote{name} and \sQuote{Euler_letter}, the latter is the notation
  from Eulers original paper; and edge attributes \code{name} (the name
  of the bridge) and \sQuote{Euler_letter}, again, Euler's notation
  from his paper.

  This dataset is in the public domain.
}
\source{Wikipedia,
  \url{http://de.wikipedia.org/wiki/K\%C3\%B6nigsberger_Br\%C3\%BCcken}
}
% \references{}
\keyword{datasets}
