\name{gclust.app}
\alias{gclust.app}
\concept{Graph Clustering}
\title{Graph Clustering Using NMF (and no SVT) -- Apparent Clusters}
\description{Performs Clustering of Graphs using Non-negative Factorization.}

\usage{
    gclust.app(glist, r=1,nmfout=FALSE,maxit=10000, nmfmethod='lee')
}

\arguments{
    \item{glist}{List of \code{igraph} Objects}
    \item{r}{Maximum Number of Clusters Allowed}
    \item{nmfout}{T/F indicating if the output from nmf should be returned; Default is FALSE}
    \item{maxit}{A number passed to irlba or svd function limiting the number of iteration; Default is 10000}
    \item{nmfmethod}{A number passed to irlba or svd, limiting the number of iterations; Default is 'lee'}
}

\value{
    \item{nmf}{An NMF object}
    \item{W}{Basis Graphs} 
    \item{H}{Probability Vector for Cluster Weights timeseries of graphs}
    \item{Xorigin}{Input Data in the matrix form}
    }

\author{Nam Lee \email{nhlee@jhu.edu}}

% \examples{
% }
