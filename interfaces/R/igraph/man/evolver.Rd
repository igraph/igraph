\name{revolver}
\alias{evolver}
\alias{revolver}
\alias{evolver.d}
\alias{revolver.d}
\alias{revolver.ad}
\alias{revolver.ade}
\alias{revolver.adi}
\alias{revolver.air}
\alias{revolver.ar}
\alias{revolver.de}
\alias{revolver.di}
\alias{revolver.dl}
\alias{revolver.e}
\alias{revolver.el}
\alias{revolver.il}
\alias{revolver.ir}
\alias{revolver.l}
\alias{revolver.r}
\alias{revolver.error.d}
\alias{revolver.error.ad}
\alias{revolver.error.ade}
\alias{revolver.error.adi}
\alias{revolver.error.air}
\alias{revolver.error.ar}
\alias{revolver.error.de}
\alias{revolver.error.di}
\alias{revolver.error.dl}
\alias{revolver.error.e}
\alias{revolver.error.el}
\alias{revolver.error.il}
\alias{revolver.error.ir}
\alias{revolver.error.l}
\alias{revolver.error.r}
\alias{revolver.d.d}
\alias{revolver.p.p}
\title{Measuring the driving force in evolving networks}
\description{These functions assume a simple evolving network model and
  measure the functional form of a so-called \emph{attractiveness
    function} governing the evolution of the network.
}
\usage{
evolver.d (nodes, kernel, outseq = NULL, outdist = NULL, m = 1, 
           directed = TRUE, verbose = igraph.par("verbose"))

revolver.d (graph, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
           verbose=igraph.par("verbose"))
revolver.ad (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(nc=2, nr=0), verbose=igraph.par("verbose"))
revolver.ade (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
                        verbose=igraph.par("verbose"))
revolver.e (graph, cats, niter=5,
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.de (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.l (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.dl (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.el (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.r (graph, window, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
           verbose=igraph.par("verbose"))
revolver.ar (graph, window, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(nc=2, nr=0), verbose=igraph.par("verbose"))
revolver.di (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.adi (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
             verbose=igraph.par("verbose"))
revolver.il (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.ir (graph, cats, window, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=igraph.par("verbose"))
revolver.air (graph, cats, window,
             niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
             verbose=igraph.par("verbose"))
revolver.d.d (graph, vtime = V(graph)$time, etime = E(graph)$time, niter = 5, 
             sd = FALSE, norm = FALSE, cites = FALSE, expected = FALSE, 
             error = TRUE, debug = matrix(nc = 2, nr = 0),
             verbose = igraph.par("verbose")) 
revolver.p.p (graph, events = get.graph.attribute(graph, "events"), 
             vtime = V(graph)$time, etime = E(graph)$time, niter = 5, sd = FALSE, 
             norm = FALSE, cites = FALSE, expected = FALSE, error = TRUE, 
             debug = matrix(nc = 2, nr = 0), verbose = igraph.par("verbose")) 
revolver.error.d (graph, kernel)
revolver.error.ad (graph, kernel)
revolver.error.ade (graph, kernel, cats)
revolver.error.adi (graph, kernel, cats)
revolver.error.air (graph, kernel, cats, window)
revolver.error.ar (graph, kernel, window)
revolver.error.de (graph, kernel, cats)
revolver.error.di (graph, kernel, cats)
revolver.error.dl (graph, kernel)
revolver.error.e (graph, kernel, cats)
revolver.error.el (graph, kernel, cats)
revolver.error.il (graph, kernel, cats)
revolver.error.ir (graph, kernel, cats, window)
revolver.error.l (graph, kernel)
revolver.error.r (graph, kernel, window)

}
\arguments{
  \item{nodes}{}
  \item{kernel}{}
  \item{outseq}{}
  \item{outdist}{}
  \item{m}{}
  \item{directed}{}
  \item{verbose}{}
  \item{graph}{}
  \item{niter}{}
  \item{sd}{}
  \item{norm}{}
  \item{cites}{}
  \item{expected}{}
  \item{error}{}
  \item{debug}{}
  \item{agebins}{}
  \item{cats}{}
  \item{window}{}
  \item{vtime}{}
  \item{etime}{}
  \item{events}{}  
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
