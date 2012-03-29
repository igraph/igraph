\name{triad.census}
\alias{triad.census}
\concept{Triad census}
\title{Triad census, subgraphs with three vertices}
\description{This function counts the different subgraphs of three
  vertices in a graph.}
\usage{
triad.census(graph)
}
\arguments{
  \item{graph}{The input graph, it should be directed. An undirected
    graph results a warning, and undefined results.}
}
\details{
  Triad census was defined by David and Leinhardt (see References
  below). Every triple of vertices (A, B, C) are classified into the 16
  possible states:
  \itemize{
    \item{003}{A,B,C, the empty graph.}
    \item{012}{A->B, C, the graph with a single directed edge.}
    \item{102}{A<->B, C, the graph with a mutual connection between two
      vertices.}
    \item{021D}{A<-B->C, the out-star.}
    \item{021U}{A->B<-C, the in-star.}
    \item{021C}{A->B->C, directed line.}
    \item{111D}{A<->B<-C.}
    \item{111U}{A<->B->C.}
    \item{030T}{A->B<-C, A->C.}
    \item{030C}{A<-B<-C, A->C.}
    \item{201}{A<->B<->C.}
    \item{120D}{A<-B->C, A<->C.}
    \item{120U}{A->B<-C, A<->C.}
    \item{120C}{A->B->C, A<->C.}
    \item{210}{A->B<->C, A<->C.}
    \item{300}{A<->B<->C, A<->C, the complete graph.}
  }

  This functions uses the RANDESU motif finder algorithm to find and
  count the subgraphs, see \code{\link{graph.motifs}}. 
}
\value{A numeric vector, the subgraph counts, in the order given in the
  above description.}
\references{
  See also Davis, J.A. and Leinhardt, S.  (1972).  The Structure of
  Positive Interpersonal Relations in Small Groups.  In J. Berger
  (Ed.), Sociological Theories in Progress, Volume 2, 218-251. 
  Boston: Houghton Mifflin.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{dyad.census}} for classifying binary relationships,
  \code{\link{graph.motifs}} for the underlying implementation.}
\examples{
g <- erdos.renyi.game(15, 45, type="gnm", directed=TRUE)
triad.census(g)
}
\keyword{graphs}
