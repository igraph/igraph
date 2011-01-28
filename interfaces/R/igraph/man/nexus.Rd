\name{nexus}
\alias{nexus}
\alias{nexus.list}
\alias{nexus.info}
\alias{nexus.get}
\alias{nexusDatasetInfo}
\alias{print.nexusDatasetInfo}
\alias{print.nexusDatasetInfoList}
\alias{summary.nexusDatasetInfoList}
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
nexus.get(id, offset=0, order=c("date", "name", "popularity"),
          nexus.url=getIgraphOpt("nexus.url"))
nexus.search(q, offset=10, limit=10,
             order=c("date", "name", "popularity"),
             nexus.url=getIgraphOpt("nexus.url"))

\method{print}{nexusDatasetInfo}(x, \dots)
\method{summary}{nexusDatasetInfoList}(object, \dots)
\method{print}{nexusDatasetInfoList}(x, \dots)
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
  \item{id}{The numeric or character id of the data set to query or
    download. Instead of the data set ids, it is possible to supply a
    \code{nexusDatasetInfo} or \code{nexusDatasetInfoList} object
    here directly and then the query is done on the corresponding data
    set(s).}
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
    \item{id}{The numeric id of the dataset.}
    \item{sid}{The character id of the dataset.}
    \item{name}{Character scalar, the name of the dataset.}
    \item{vertices/edges}{Character, the number of vertices and edges in
      the graph(s). Vertices and edges are separated by a slash, and
      if the data set consists of multiple networks, then they are
      separated by spaces.}
    \item{tags}{Character vector, the tags of the dataset. Directed
      graph have the tags \sQuote{directed}. Undirected graphs are
      tagged as \sQuote{undirected}. Other common tags are:
      \sQuote{weighted}, \sQuote{bipartite}, \sQuote{social network},
      etc.}
    \item{networks}{The ids and names of the networks in the data
      set. The numeric and character id are separated by a slash, and
      multiple networks are separated by spaces.}
  }

  \code{nexusDatasetInfo} objects returned by \code{nexus.info} have the
  following additional fields:
  \describe{
    \item{date}{Character scalar, e.g. \sQuote{2011-01-09}, the date
      when the dataset was added to the database.}
    \item{formats}{Character vector, the data formats in which the data
      set is available. The various formats are separated by semicolons.}
    \item{licence}{Character scalar, the licence of the dataset.}
    \item{licence url}{Character scalar, the URL of the licence of the
      dataset. Pleaase make sure you consult this before using a
      dataset.}
    \item{summary}{Character scalar, the short description of
      the dataset, this is usually a single sentence.}
    \item{description}{Character scalar, the full description of the
      dataset.}
    \item{citation}{Character scalar, the paper(s) describing the
      dataset. Please cite these papers if you are using the dataset in
      your research, the licence of most datasets requires this.}
    \item{attributes}{A list of lists, each list entry is a graph,
      vertex or edge attribute and has the following entries:
      \describe{
	\item{type}{Type of the attribute, either \sQuote{graph},
	  \sQuote{vertex} or \sQuote{edge}.}
	\item{datatype}{Data type of the attribute, currently it can be
	  \sQuote{numeric} and \sQuote{string}.}
	\item{name}{Character scalar, the name of the attribute.}
	\item{description}{Character scalar, the description of the
	  attribute.}
      }
    }
  }

  The results of the Nexus queries are printed to the screen in a
  consise format, similar to the format of igraph graphs. A data set
  list (typically the result of \code{nexus.list} and
  \code{nexus.search}) looks like this: \preformatted{NEXUS 1-5/18 -- data set list
[1] kaptail.4         39/109-223   #18 Kapferer tailor shop
[2] condmatcollab2003 31163/120029 #17 Condensed matter collaborations+
[3] condmatcollab     16726/47594  #16 Condensed matter collaborations+
[4] powergrid         4941/6594    #15 Western US power grid
[5] celegansneural    297/2359     #14 C. Elegans neural network }
  Each line here represents a data set, and the following information is
  given about them: the character id of the data set
  (e.g. \code{kaptail} or \code{powergrid}), the number of vertices and
  number of edges in the graph of the data sets. For data sets with
  multiple graphs, intervals are given here. Then the numeric id of the
  data set and the reamining space is filled with the name of the data
  set.

  Summary information about an individual Nexus data set is printed
  as \preformatted{NEXUS B--- 39 109-223 -- #18 Kapferer tailor shop
+ tags: directed; social network; undirected
+ nets: 1/KAPFTI2; 2/KAPFTS2; 3/KAPFTI1; 4/KAPFTS1}
  This is very similar to the header that is used for printing igraph
  graphs, but there are some differences as well. The four characters
  after the \code{NEXUS} word give the most important properties of the
  graph(s): the first is \sQuote{\code{U}} for undirected and
  \sQuote{\code{D}} for directed graphs, and \sQuote{\code{B}} if the
  data set contains both directed and undirected graphs. The second is
  \sQuote{\code{N}} named graphs. The third character is
  \sQuote{\code{W}} for weighted graphs, the fourth is \sQuote{\code{B}}
  if the data set contains bipartite graphs. The end of the first list
  contains the numeric id and the name of the data set. The second row
  lists the data set tags, and the third row the networks that are
  included in the data set.

  Detailed data set information is printed similarly, but it contains
  more fields.
}
\value{
  \code{nexus.list} and \code{nexus.search} return a list of
  \code{nexusDatasetInfo} objects. The list also has these attributes:
  \describe{
    \item{size}{The number of data sets returned by the query.}
    \item{totalsize}{The total number of data sets found for the query.}
    \item{offset}{The offset parameter of the query.}
    \item{limit}{The limit parameter of the query.}
  }
  
  \code{nexus.info} returns a single \code{nexusDatasetInfo} object.

  \code{nexus.get} returns an igraph graph object, or a list of graph
  objects, if the data set consists of multiple networks.
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

## Any word can match
nexus.search("blog or US or karate")
}
}
\keyword{graphs}
