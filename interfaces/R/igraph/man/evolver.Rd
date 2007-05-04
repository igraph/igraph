\name{evolver}
\alias{evolver}
\alias{evolver.d}
\alias{evolver.ad}
\alias{evolver.ade}
\alias{evolver.adi}
\alias{evolver.air}
\alias{evolver.ar}
\alias{evolver.de}
\alias{evolver.di}
\alias{evolver.dl}
\alias{evolver.e}
\alias{evolver.el}
\alias{evolver.il}
\alias{evolver.ir}
\alias{evolver.l}
\alias{evolver.r}
\title{Measuring the driving force in evolving networks}
\description{These functions assume a simple evolving network model and
  measure the functional form of a so-called \emph{attractiveness
    function} governing the evolution of the network.
}
\usage{
evolver.d (graph, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
           verbose=igraph.par("verbose"))
evolver.ad (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(nc=2, nr=0), verbose=igraph.par("verbose"))
evolver.ade (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
                        verbose=igraph.par("verbose"))
evolver.e (graph, cats, niter=5,
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.de (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.l (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.dl (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.el (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.r (graph, window, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
           verbose=igraph.par("verbose"))
evolver.ar (graph, window, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(nc=2, nr=0), verbose=igraph.par("verbose"))
evolver.di (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.adi (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
             verbose=igraph.par("verbose"))
evolver.il (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.ir (graph, cats, window, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
evolver.air (graph, cats, window,
             niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
             verbose=igraph.par("verbose"))
}
\arguments{
  \item{graph}{}
  \item{niter}{}
  \item{sd}{}
  \item{norm}{}
  \item{cites}{}
  \item{expected}{}
  \item{error}{}
  \item{debug}{}
  \item{verbose}{}
  \item{agebins}{}
  \item{cats}{}
  \item{window}{}
}
\details{
  The functions should be considered as experimental, so no
  detailed documentation yet. Sorry.
}
\value{
  A named list.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% seealso{}
% \examples{}
\keyword{graphs}
