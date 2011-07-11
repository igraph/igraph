\name{USairports}
\docType{data}
\alias{USairports}
\title{US airport network, 2010 December}
\description{
  The network of passanger flights between airports in the United
  States. The data set was compiled based on flights in 2010
  December. This network is directed and edge directions correspond to
  flight directions. Each edge is specific to a single carrier aircraft
  type. Multiple carriers between the same two airports are denoted by
  multiple edges.

  See information about the included meta-data below.
}
\usage{
USairports
}
\format{
  A directed \code{igraph} graph object, with multiple edges. It has a
  \sQuote{name} graph attribute, and several vertex and edge
  attributes. The vertex attributes:
  \describe{
    \item{name}{Symbolic vertex name, this is the three letter IATA
      airport code.}
    \item{City}{City and state, where the airport is located.}
    \item{Position}{Position of the airport, in WGS coordinates.}
  }

  Edge attributes:
  \describe{
    \item{Carrier}{Name of the airline. The network includes both
      domestic and international carriers that performed at least one
      flight in December of 2010.}
    \item{Departures}{The number of departures (for a given airline and
      aircraft type.}
    \item{Seats}{The total number of seats available on the flights
      carried out by a given airline, using a given aircraft type.}
    \item{Passengers}{The total number of passangers on the flights
      carried out by a given airline, using a given aircraft type.}
    \item{Aircraft}{Type of the aircraft.}
    \item{Distance}{The distance between the two airports, in miles.}
  }
}
\source{
  Most of this information was downloaded from The Research and
  Innovative Technology Administration (RITA). See
  \url{http://www.rita.dot.gov/about_rita/} for details. The airport
  position information was collected from Wikipedia and other public
  online sources.
}
% \references{}
\keyword{datasets}
