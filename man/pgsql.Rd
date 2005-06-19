\name{pgsql}
\alias{pgsql}
\alias{pgsql.list.graphs}
\alias{pgsql.get.graph}
\alias{pgsql.clear}
\title{PGSQL specific functions}
\description{}
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
  \item{host}{}
  \item{dbname}{}
  \item{user}{}
  \item{gname}{}
  \item{all}{}
}
\details{}
\value{}
\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{}
\examples{}
\keyword{graphs}

  