\name{UKfaculty}
\docType{data}
\alias{UKfaculty}
\title{Friendship network of a UK university faculty}
\description{The personal friendship network of a faculty of a UK
  university, consisting of 81 vertices (individuals) and 817 directed
  and weighted connections. The school affiliation of each individual is
  stored as a vertex attribute. This dataset can serve as a testbed for
  community detection algorithms.}
\usage{
UKfaculty
}
\format{
  A directed \code{igraph} graph object with vertex attribute
  \sQuote{Group}, the numeric id of the school affiliation, and edge
  attribute \sQuote{weight}, i.e. the graph is weighted.

  This dataset is licensed under a Creative Commons
  Attribution-Share Alike 2.0 UK: England & Wales License, 
  see \url{http://creativecommons.org/licenses/by-sa/2.0/uk/} for details.
  Please cite the reference below if you use this dataset.  
}
\source{See reference below.}
\references{Nepusz T., Petroczi A., Negyessy L., Bazso F.: Fuzzy
  communities and the concept of bridgeness in complex
  networks. Physical Review E 77:016107, 2008.
}
\keyword{datasets}
