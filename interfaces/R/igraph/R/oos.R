## Performs out-of-sample embedding for a graph given 
## some configuration for the in-sample points

## The input is a matrix A12 of size nxm where
## n = # of in-sample points
## m = # of out-of-sample points
## and a matrix Xhat of size nxd containing
## the embedding of the in-sample points.
## The output is a matrix of size mxd containing the
## embedding of the out-of-sample points.

## The oos algorithm itself is nothing more than a least square regression, i.e., 
## Yhat = \argmin \| A12 - Xhat Y^{\top} \|_{F}
oos <- function(A12, Xhat){
    Yhat <- t(Matrix::solve(t(Xhat)%*%Xhat, t(Xhat) %*% A12))
    return(Yhat)
}
