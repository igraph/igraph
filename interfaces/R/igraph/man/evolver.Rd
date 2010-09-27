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
\alias{revolver.ml.AD.alpha.a.beta}
\alias{revolver.ml.AD.dpareto}
\alias{revolver.ml.AD.dpareto.eval}
\alias{revolver.ml.ad}
\alias{revolver.ml.ade}
\alias{revolver.ml.ADE.alpha.a.beta}
\alias{revolver.ml.ADE.dpareto}
\alias{revolver.ml.ADE.dpareto.eval}
\alias{revolver.ml.ADE.dpareto.evalf}
\alias{revolver.ml.d}
\alias{revolver.ml.D.alpha}
\alias{revolver.ml.D.alpha.a}
\alias{revolver.ml.de}
\alias{revolver.ml.DE.alpha.a}
\alias{revolver.ml.df}
\alias{revolver.ml.f}
\alias{revolver.ml.l}
\alias{revolver.probs.ad}
\alias{revolver.probs.ade}
\alias{revolver.probs.ADE.dpareto}
\alias{revolver.probs.d}
\alias{revolver.probs.de}
\title{Measuring the driving force in evolving networks}
\description{These functions assume a simple evolving network model and
  measure the functional form of a so-called \emph{attractiveness
    function} governing the evolution of the network.
}
\usage{
evolver.d (nodes, kernel, outseq = NULL, outdist = NULL, m = 1, 
           directed = TRUE, verbose = getIgraphOpt("verbose"))

revolver.d (graph, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
           verbose=getIgraphOpt("verbose"))
revolver.ad (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(nc=2, nr=0), verbose=getIgraphOpt("verbose"))
revolver.ade (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
                        verbose=getIgraphOpt("verbose"))
revolver.e (graph, cats, niter=5, st=FALSE,
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.de (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.l (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.dl (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.el (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.r (graph, window, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric(),
           verbose=getIgraphOpt("verbose"))
revolver.ar (graph, window, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(nc=2, nr=0), verbose=getIgraphOpt("verbose"))
revolver.di (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.adi (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
             verbose=getIgraphOpt("verbose"))
revolver.il (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.ir (graph, cats, window, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric(), verbose=getIgraphOpt("verbose"))
revolver.air (graph, cats, window,
             niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(nc=2, nr=0),
             verbose=getIgraphOpt("verbose"))
revolver.d.d (graph, vtime = V(graph)$time, etime = E(graph)$time, niter = 5, 
             sd = FALSE, norm = FALSE, cites = FALSE, expected = FALSE, 
             error = TRUE, debug = matrix(nc = 2, nr = 0),
             verbose = getIgraphOpt("verbose")) 
revolver.p.p (graph, events = get.graph.attribute(graph, "events"), 
             vtime = V(graph)$time, etime = E(graph)$time, niter = 5, sd = FALSE, 
             norm = FALSE, cites = FALSE, expected = FALSE, error = TRUE, 
             debug = matrix(nc = 2, nr = 0), verbose = getIgraphOpt("verbose")) 
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

revolver.ml.ade (graph, niter, cats, agebins = 300, delta = 1e-10,
    filter = NULL, verbose = getIgraphOpt("verbose"))
revolver.ml.d (graph, niter, delta = 1e-10, filter = NULL,
    verbose = getIgraphOpt("verbose")) 
revolver.ml.de (graph, niter, cats, delta = 1e-10, filter = NULL,
    verbose = getIgraphOpt("verbose")) 
revolver.ml.df (graph, niter, delta = 1e-10,
    verbose = getIgraphOpt("verbose")) 
revolver.ml.f (graph, niter, delta = 1e-10,
    verbose = getIgraphOpt("verbose")) 
revolver.ml.l (graph, niter, agebins = 300, delta = 1e-10)

revolver.ml.AD.alpha.a.beta (graph, alpha, a, beta, abstol = 1e-08,
    reltol = 1e-08, maxit = 1000, agebins = 300, filter = NULL,
    verbose = getIgraphOpt("verbose"))
revolver.ml.AD.dpareto (graph, alpha, a, paralpha, parbeta, parscale,
    abstol = 1e-08, reltol = 1e-08, maxit = 1000, agebins = 300, filter
    = NULL, verbose = getIgraphOpt("verbose"))
revolver.ml.ADE.alpha.a.beta (graph, cats, alpha, a, beta, coeffs,
    abstol = 1e-08, reltol = 1e-08, maxit = 1000, agebins = 300, filter
    = NULL)
revolver.ml.ADE.dpareto (graph, cats, alpha, a, paralpha, parbeta, parscale, 
    coeffs, abstol = 1e-08, reltol = 1e-08, maxit = 1000, agebins = 300, 
    filter = NULL)
revolver.ml.D.alpha (graph, alpha, abstol = 1e-08, reltol = 1e-08, maxit
    = 1000, filter = NULL, verbose = getIgraphOpt("verbose"))
revolver.ml.D.alpha.a (graph, alpha, a, abstol = 1e-08, reltol = 1e-08,
    maxit = 1000, filter = NULL, verbose = getIgraphOpt("verbose"))
revolver.ml.DE.alpha.a (graph, cats, alpha, a, coeffs, abstol = 1e-08,
    reltol = 1e-08, maxit = 1000, filter = NULL)

revolver.ml.AD.dpareto.eval (graph, alpha, a, paralpha, parbeta,
    parscale, agebins = 300, filter = NULL) 
revolver.ml.ADE.dpareto.eval (graph, cats, alpha, a, paralpha, parbeta,
    parscale, coeffs, agebins = 300, filter = NULL) 
revolver.ml.ADE.dpareto.evalf (graph, cats, par, agebins, filter = NULL)

revolver.probs.ad (graph, kernel, ntk = FALSE)
revolver.probs.ade (graph, kernel, cats)
revolver.probs.d (graph, kernel, ntk = FALSE)
revolver.probs.de (graph, kernel, cats)
revolver.probs.ADE.dpareto (graph, par, cats, gcats, agebins)
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
  \item{delta}{}
  \item{filter}{}
  \item{alpha}{}
  \item{a}{}
  \item{paralpha}{}
  \item{parbeta}{}
  \item{parscale}{}
  \item{abstol}{}
  \item{reltol}{}
  \item{maxit}{}
  \item{beta}{}
  \item{coeffs}{}
  \item{par}{}
  \item{ntk}{}
  \item{gcats}{}
  \item{st}{}
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
