\name{pgsql}
\alias{pgsql}
\alias{pgsql.list.graphs}
\alias{pgsql.get.graph}
\alias{pgsql.clear}
\title{PGSQL specific functions}
\description{The \code{plsql} storage type is a toy implementation of a
  database-based storage type where the graphs are stored in an SQL
  database. It is too slow for production use. This manual page
  documents some \code{pgsql} storage type specific functions.}
\usage{
pgsql.list.graphs(host = igraph.par("pgsql.host"),
                  dbname = igraph.par("pgsql.dbname"),
                  user = igraph.par("pgsql.user"))
pgsql.get.graph(gname, host = igraph.par("pgsql.host"),
                dbname = igraph.par("pgsql.dbname"),
                user = igraph.par("pgsql.user"), \dots)
pgsql.clear(all = FALSE, host = igraph.par("pgsql.host"),
            dbname = igraph.par("pgsql.dbname"), 
            user = igraph.par("pgsql.user"))
}
\arguments{
  \item{host}{The host on which the database resides.}
  \item{dbname}{The name of the database.}
  \item{user}{The name of the user.}
  \item{gname}{Character constant, the name of the graph.}
  \item{all}{Logical constant, whether to erase all graphs from the
    database.}
}
\details{
  \code{pgsql.list.graphs} simply lists the ids of the graphs in a
  database.

  \code{pgsql.get.graph} gives back the graph with the given id from
  a database.

  \code{pgsql.clear} erases a single graph, or optionally all graphs
  from a database.
}
\value{
  For \code{pgsql.list.graphs} a vector with the graph ids in the
  database.

  For \code{pgsql.get.graph} a graph object, with storage type
  \code{pgsql}.

  For \code{pgsql.clear returns nothing.}
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
\dontrun{g <- graph.ring(10, type="pgsql", host="localhost",
                         dbname="igraph", user="igraph")}
\dontrun{pgsql.list.graphs(host="localhost", dbname="igraph", user="igraph")}
}
\keyword{graphs}
