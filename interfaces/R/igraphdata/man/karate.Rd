\name{karate}
\docType{data}
\alias{karate}
\title{Zachary's karate club network}
\description{
  Social network between members of a university karate club, led by
  president John A. and karate instructor Mr. Hi (pseudonyms).

  The edge weights are the number of common activities the club
  members took part of. These activities were:
  \enumerate{
    \item Association in and between academic classes at the university.
    \item Membership in Mr. Hi's private karate studio on the east side
    of the city where Mr. Hi taught nights as a part-time instructor.
    \item Membership in Mr. Hi's private karate studio on the east side
    of the city, where many of his supporters worked out on weekends. 
    \item Student teaching at the east-side karate studio referred to in
    (2). This is different from (2) in that student teachers interacted
    with each other, but were prohibited from interacting with their
    students. 
    \item Interaction at the university rathskeller, located in the same
    basement as the karate club's workout area.
    \item Interaction at a student-oriented bar located across the
    street from the university campus.
    \item Attendance at open karate tournaments held through the area at
    private karate studios.
    \item Attendance at intercollegiate karate tournaments held at local
    universities. Since both open and intercollegiate tournaments were
    held on Saturdays, attendance at both was impossible.
  }

  Zachary studied conflict and fission in this network, as the karate
  club was split into two separate clubs, after long disputes between
  two factions of the club, one led by John A., the other by Mr. Hi.

  The \sQuote{Faction} vertex attribute gives the faction memberships of
  the actors. After the split of the club, club members chose their
  new clubs based on their factions, except actor no. 9, who was in John
  A.'s faction but chose Mr. Hi's club.
}
\usage{
karate
}
\format{
  An undirected \code{igraph} graph object. Vertex no. 1 is Mr. Hi,
  vertex no. 34 corresponds to John A.

  Graph attributes: \sQuote{name}, \sQuote{Citation}, \sQuote{Author}.

  Vertex attribute: \sQuote{Faction}.

  Edge attribute: \sQuote{weight}.
}
\source{
  See reference below.
}
\references{
  Wayne W. Zachary. An Information Flow Model for Conflict and Fission
  in Small Groups. \emph{Journal of Anthropological Research} Vol. 33,
  No. 4 452-473
}
\keyword{datasets}
