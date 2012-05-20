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
           directed = TRUE)

revolver.d (graph, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric())
revolver.ad (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(ncol=2, nrow=0))
revolver.ade (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(ncol=2, nrow=0))
revolver.e (graph, cats, niter=5, st=FALSE,
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric())
revolver.de (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric())
revolver.l (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
           sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
           error=TRUE, debug=numeric())
revolver.dl (graph, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric())
revolver.el (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric())
revolver.r (graph, window, niter=5, sd=FALSE, norm=FALSE,
           cites=FALSE, expected=FALSE, error=TRUE, debug=numeric())
revolver.ar (graph, window, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE, error=TRUE,
            debug=matrix(ncol=2, nrow=0))
revolver.di (graph, cats, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric())
revolver.adi (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(ncol=2, nrow=0))
revolver.il (graph, cats, niter=5, agebins=max(vcount(graph)/7100, 10),
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric())
revolver.ir (graph, cats, window, niter=5,
            sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
            error=TRUE, debug=numeric())
revolver.air (graph, cats, window,
             niter=5, agebins=max(vcount(graph)/7100, 10),
             sd=FALSE, norm=FALSE, cites=FALSE, expected=FALSE,
             error=TRUE, debug=matrix(ncol=2, nrow=0))
revolver.d.d (graph, vtime = V(graph)$time, etime = E(graph)$time, niter = 5, 
             sd = FALSE, norm = FALSE, cites = FALSE, expected = FALSE, 
             error = TRUE, debug = matrix(ncol = 2, nrow = 0))
revolver.p.p (graph, events = get.graph.attribute(graph, "events"), 
             vtime = V(graph)$time, etime = E(graph)$time, niter = 5, sd = FALSE, 
             norm = FALSE, cites = FALSE, expected = FALSE, error = TRUE, 
             debug = matrix(ncol = 2, nrow = 0)) 
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
    filter = NULL)
revolver.ml.d (graph, niter, delta = 1e-10, filter = NULL)
revolver.ml.de (graph, niter, cats, delta = 1e-10, filter = NULL)
revolver.ml.df (graph, niter, delta = 1e-10)
revolver.ml.f (graph, niter, delta = 1e-10)
revolver.ml.l (graph, niter, agebins = 300, delta = 1e-10)

revolver.ml.AD.alpha.a.beta (graph, alpha, a, beta, abstol = 1e-08,
    reltol = 1e-08, maxit = 1000, agebins = 300, filter = NULL)
revolver.ml.AD.dpareto (graph, alpha, a, paralpha, parbeta, parscale,
    abstol = 1e-08, reltol = 1e-08, maxit = 1000, agebins = 300, filter
    = NULL)
revolver.ml.ADE.alpha.a.beta (graph, cats, alpha, a, beta, coeffs,
    abstol = 1e-08, reltol = 1e-08, maxit = 1000, agebins = 300, filter
    = NULL)
revolver.ml.ADE.dpareto (graph, cats, alpha, a, paralpha, parbeta, parscale, 
    coeffs, abstol = 1e-08, reltol = 1e-08, maxit = 1000, agebins = 300, 
    filter = NULL)
revolver.ml.D.alpha (graph, alpha, abstol = 1e-08, reltol = 1e-08, maxit
    = 1000, filter = NULL)
revolver.ml.D.alpha.a (graph, alpha, a, abstol = 1e-08, reltol = 1e-08,
    maxit = 1000, filter = NULL)
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
  \item{nodes}{The number of vertices in the generated network.}
  \item{kernel}{The kernel function, a vector, matrix or array,
    depending on the number of model parameters.}
  \item{outseq}{The out-degree sequence, or \code{NULL} if no out-degree
    sequence is used.}
  \item{outdist}{The out-degree distribution, or \code{NULL} if all
    vertices have the same out-degree. This argument is ignored if the
    \code{outseq} argument is not \code{NULL}.}
  \item{m}{Numeric scalar, the out-degree of the verticec. It is ignored
    if at least one of \code{outseq} and \code{outdist} is not
    \code{NULL}.}
  \item{directed}{Logical scalar, whether to create a directed graph.}
  \item{graph}{The input graph.}
  \item{niter}{The number of iterations to perform.}
  \item{sd}{Logical scalar, whether to return the standard deviation of
    the estimates.}
  \item{norm}{Logical scalar, whether to return the normalizing factors.}
  \item{cites}{Logical scalar, whether to return the number of citations
    to the different vertex types.}
  \item{expected}{Logical scalar, whether to return the expected number
    of citations for the different vertex types.}
  \item{error}{Logical scalar, whether to return the error of the fit.}
  \item{debug}{Currently not used.}
  \item{agebins}{The number of bins for vertex age.}
  \item{cats}{The number of categories to use.}
  \item{window}{The width of the time window to use, measured in number
    of vertices.}
  \item{vtime}{Numeric vector, the time steps when the vertices where
    added to the network.}
  \item{etime}{Numeric vector, the time steps when the edges where added
    to the network.}
  \item{events}{A list of numeric vectors, each vector represents an
    event, with the participation of the listed vertices.}
  \item{delta}{Real scalar, the error margin that is allowed for the
    convergence.}
  \item{filter}{Logical vector, length is the number of vertices. Only
    vertices corresponding to \code{TRUE} entries are used in the
    fitting.}
  \item{alpha}{Starting value for the \sQuote{\code{alpha}} parameter.}
  \item{a}{Starting value for the \sQuote{\code{a}} parameter.}
  \item{paralpha}{Starting value for the \sQuote{\code{paralpha}}
    (Pareto alpha) parameter.}
  \item{parbeta}{Starting value for the \sQuote{\code{parbeta}} (Pareto
    beta) parameter.}
  \item{parscale}{Starting value for the \sQuote{\code{parscale}}
    (Pareto scale) parameter.}
  \item{abstol}{Real scalar, absolute tolerance for the ML fitting.}
  \item{reltol}{Real scalar, relative tolerance for the ML fitting.}
  \item{maxit}{Numeric scalar, the maximum number of iterations.}
  \item{beta}{Real scalar, starting value for the \sQuote{\code{beta}}
    parameter.}
  \item{coeffs}{Numeric vector, starting values for the coefficients.}
  \item{par}{Pareto parameters for the different vertex types, in a
    matrix.}
  \item{ntk}{Logical scalar, whether to return the Ntk values.}
  \item{gcats}{Numeric vector, the vertex types.}
  \item{st}{Logical scalar, whether to return the S(t) values.}
}
\details{
  The functions should be considered as experimental, so no
  detailed documentation yet. Sorry.
}
\value{
  A named list.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% seealso{}
% \examples{}
\keyword{graphs}
