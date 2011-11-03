\name{yeast}
\docType{data}
\alias{yeast}
\title{Yeast protein interaction network}
\description{
  Comprehensive protein-protein interaction maps promise to reveal many
  aspects of the complex regulatory network underlying cellular
  function.

  This data set was compiled by von Mering et al. (see reference below),
  combining various sources. Only the interactions that have 
  \sQuote{high} and \sQuote{medium} confidence are included here.
}
\usage{
yeast
}
\format{
  An undirected \code{igraph} graph object. Its graph attributes:
  \sQuote{name}, \sQuote{Citation}, \sQuote{Author},
  \sQuote{URL}. \sQuote{Classes}. The \sQuote{Classes}
  attribute contain the key for the classification labels of the
  proteins, in a data frame, the original MIPS categories are given
  after the semicolon:
  \describe{
    \item{E}{energy production; energy}
    \item{G}{aminoacid metabolism; aminoacid metabolism} 
    \item{M}{other metabolism; all remaining metabolism categories}
    \item{P}{translation; protein synthesis}
    \item{T}{transcription; transcription, but without subcategory
      \sQuote{transcriptional control}}
    \item{B}{transcriptional control; subcategory
      \sQuote{transcriptional control}}
    \item{F}{protein fate; protein fate (folding, modification,
      destination)}
    \item{O}{cellular organization; cellular transport and transport
      mechanisms}
    \item{A}{transport and sensing; categories \sQuote{transport
	facilitation} and \sQuote{regulation of / interaction with
	cellular environment}}
    \item{R}{stress and defense; cell rescue, defense and virulence}
    \item{D}{genome maintenance; DNA processing and cell cycle}
    \item{C}{cellular fate / organization; categories \sQuote{cell fate}
      and \sQuote{cellular communication /  signal transduction} and
      \sQuote{control of cellular organization}}
    \item{U}{uncharacterized; categories \sQuote{not yet clear-cut} and
      \sQuote{uncharacterized}}
  }
  
  Vertex attributes: \sQuote{name}, \sQuote{Description},
  \sQuote{Class}, the last one contains the class of the protein,
  accoring to the classification above.

  Note that some proteins in the network did not appear in the
  annotation files, the \sQuote{Class} and \sQuote{Description}
  attributes are \code{NA} for these.
}
\source{The data was downloaded from
  \url{http://www.nature.com/nature/journal/v417/n6887/suppinfo/nature750.html}.
}
\references{
Comparative assessment of large-scale data sets of protein-protein
  interactions. Christian von Mering, Roland Krause, Berend Snel,
  Michael Cornell, Stephen G. Oliver, Stanley Fields and Peer
  Bork. \emph{Nature} 417, 399-403 (2002) 
}
\keyword{datasets}
