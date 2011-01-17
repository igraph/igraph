\name{nexus}
\alias{nexus}
\alias{nexus.list}
\alias{nexus.info}
\alias{nexus.get}
\alias{nexusDatasetInfo}
\alias{print.nexusDatasetInfo}
\concept{The Nexus network repository}
\title{Query and download from the Nexus network repository}
\description{The Nexus network repository is an online collection of
  network data sets. These functions can be used to query it and
  download data from it, directly as an igraph graph.}
\usage{
nexus.list(tags=NULL, offset=0, limit=10,
           operator=c("or", "and"),
           order=c("date", "name", "popularity"),
           nexus.url=getIgraphOpt("nexus.url"))
nexus.info(id, nexus.url=getIgraphOpt("nexus.url"))
nexus.get(id, q, offset=0, order=c("date", "name", "popularity"),
          nexus.url=getIgraphOpt("nexus.url"))
nexus.search(q, offset=10, limit=10,
             order=c("date", "name", "popularity"),
             nexus.url=getIgraphOpt("nexus.url"))

\method{print}{nexusDatasetInfo}(x, \dots)
}
\arguments{
  \item{tags}{A character vector, the tags that are searched. If not
    given (or \code{NULL}), then all datasets are listed.}
  \item{offset}{An offset to select part of the results. Results are
    listed from \code{offset}+1.}
  \item{limit}{The maximum number of results to return.}
  \item{operator}{A character scalar. If \sQuote{or} (the default),
    then all datasets that have at least one of the given tags, are
    returned. If it if \sQuote{and}, then only datasets that have all
    the given tags, are returned.}
  \item{order}{The ordering of the results, possible values are:
    \sQuote{date}, \sQuote{name}, \sQuote{popularity}.}
  \item{id}{The numeric id of the data set to query or download.}
  \item{q}{Nexus search string. See examples below. For the complete
    documentation please see the Nexus homepage at
    \url{http://nexus.igraph.org}.}
  \item{nexus.url}{The URL of the Nexus server. Don't change this from
    the default, unless you set up your own Nexus server.}
  \item{x}{The \code{nexusDatasetInfo} object to print.}
  \item{\dots}{Currently ignored.}
}
\details{
  Nexus is an online repository of networks, with an API that allow
  programatic queries against it, and programatic data download as
  well.

  The \code{nexus.list} and \code{nexus.info} functions query the
  online database. They both return \code{nexusDatasetInfo} objects.
  \code{nexus.info} returns more information than \code{nexus.list}.

  \code{nexus.search} searches Nexus, and returns a list of data sets,
  as \code{nexusDatasetInfo} objects. See below for some search
  examples.

  \code{nexus.get} downloads a data set from Nexus, based on its numeric
  id, or based on a Nexus search string. For search strings, only the
  first search hit is downloaded, but see also the \code{offset}
  argument. (If there are not data sets found, then the function returns
  an error.)
  
  The \code{nexusDatasetInfo} objects returned by \code{nexus.list} have
  the following fields:
  \describe{
    \item{Id}{The numeric id of the dataset.}
    \item{Name}{Character scalar, the name of the dataset.}
    \item{Vertices}{Numeric scalar, the number of vertices in the
      graph.}
    \item{Edges}{Numeric scalar, the number of edges in the graph.}
    \item{Tags}{Character vector, the tags of the dataset. Directed
      graph have the tags \sQuote{directed}. Undirected graphs are
      tagged as \sQuote{undirected}. Other common tags are:
      \sQuote{weighted}, \sQuote{bipartite}, \sQuote{social network},
      etc.}
  }

  \code{nexusDatasetInfo} objects returned by \code{nexus.info} have the
  following additional fields:
  \describe{
    \item{Date}{Character scalar, e.g. \sQuote{2011-01-09}, the date
      when the dataset was added to the database.}
    \item{Licence}{Character scalar, the licence of the dataset.}
    \item{Licence URL}{Character scalar, the URL of the licence of the
      dataset. Pleaase make sure you consult this before using a
      dataset.}
    \item{Short description}{Character scalar, the short description of
      the dataset, this is usually a single sentence.}
    \item{Description}{Character scalar, the full description of the
      dataset.}
    \item{Citation}{Character scalar, the paper(s) describing the
      dataset. Please cite these papers if you are using the dataset in
      your research, the licence of most datasets requires this.}
    \item{Attributes}{A list of lists, each list entry is a graph,
      vertex or edge attribute and has the following entries:
      \describe{
	\item{Type}{Type of the attribute, either \sQuote{graph},
	  \sQuote{vertex} or \sQuote{edge}.}
	\item{Datatype}{Data type of the attribute, currently it can be
	  \sQuote{numeric} and \sQuote{string}.}
	\item{Name}{Character scalar, the name of the attribute.}
	\item{Description}{Character scalar, the description of the
	  attribute.}
      }
    }
  }
}
\value{
  \code{nexus.list} and \code{nexus.search} return a list of
  \code{nexusDatasetInfo} objects. The list also has these attributes:
  \describe{
    \item{Size}{The number of data sets returned by the query.}
    \item{Totalsize}{The total number of data sets found for the query.}
    \item{Offset}{The offset parameter of the query.}
    \item{Limit}{The limit parameter of the query.}
  }
  
  \code{nexus.info} returns a single \code{nexusDatasetInfo} object.

  \code{nexus.get} returns an igraph graph object.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{ \url{http://nexus.igraph.org} }
\examples{
\dontrun{nexus.list(tag="weighted")
nexus.list(limit=3, order="name")
nexus.info(2)
g <- nexus.get(2)
summary(g)

## Data sets related to 'US':
nexus.search("US")

## Search for data sets that have 'network' in their name:
nexus.search("name:network")

## Download the first data set that matches
nexus.get(q="Euler")

## Any word can match
nexus.search("blog or US or karate")
}
}
\keyword{graphs}
