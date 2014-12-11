
solve_LSAP <- function (x, maximum = FALSE) {
  if (!is.matrix(x) || any(x < 0)) {
    stop("x must be a matrix with nonnegative entries.")
  }
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr > nc) stop("x must not have more rows than columns.")
  if (nc > nr)  x <- rbind(x, matrix(2 * sum(x), nc - nr, nc))
  if (maximum)  x <- max(x) - x
  storage.mode(x) <- "double"
  out <- .Call("R_igraph_solve_lsap", x, as.integer(nc),
               PACKAGE = "igraph") + 1L
  out[seq_len(nr)]
}

match_vertices <- function(A, B, m, start, iteration) {
  ## Seeds are assumed to be vertices 1:m in both graphs
  totv <- ncol(A)
  n <- totv - m
  if (m != 0) {
    A12 <- A[1:m, (m+1):(m+n), drop=FALSE]
    A21 <- A[(m+1):(m+n), 1:m, drop=FALSE]
    B12 <- B[1:m, (m+1):(m+n), drop=FALSE]
    B21 <- B[(m+1):(m+n), 1:m, drop=FALSE]
  }
  if ( m==0 ) {
    A12 <- Matrix::Matrix(0, n, n)
    A21 <- Matrix::Matrix(0, n, n)
    B12 <- Matrix::Matrix(0, n, n)
    B21 <- Matrix::Matrix(0, n, n)
  }
  A22 <- A[(m+1):(m+n), (m+1):(m+n)]
  B22 <- B[(m+1):(m+n), (m+1):(m+n)]
  patience <- iteration
  tol <- 1
  P <- start
  toggle <- 1
  iter <- 0
  while (toggle == 1 & iter < patience)  {
    iter <- iter+1
    x <-  A21 %*% Matrix::t(B21)
    y <-  Matrix::t(A12) %*% B12
    z <-  A22 %*% P %*% Matrix::t(B22)
    w <-  Matrix::t(A22) %*% P %*% B22
    Grad <- x + y + z + w
    ind <- unclass(solve_LSAP(as.matrix(Grad), maximum = TRUE))
    ind2 <- cbind(1:n, ind)
    T <- Matrix::Diagonal(n)
    T <- T[ind, ]
    wt <- Matrix::t(A22)[,order(ind)] %*% B22
    c <- sum(w * P)
    d <- sum(wt * P) + sum(w [ ind2 ])
    e <- sum(wt[ind2])
    u <- sum(P * (x + y))
    v <- sum((x + y)[ind2])
    if ( c-d+e == 0 && d-2*e+u-v == 0) {
      alpha <- 0
    } else {
      alpha <- -(d-2*e+u-v) / (2*(c-d+e))}
    f0 <- 0
    f1 <-  c-e+u-v
    falpha <- (c-d+e) * alpha^2 + (d-2*e+u-v) * alpha
    if (alpha < tol && alpha > 0 && falpha > f0 && falpha > f1) {
      P <-  alpha*P + (1-alpha) * T
    } else if (f0 > f1) {
      P <- T
    } else {
      toggle <- 0
    }
  }
  D <- P
  corr <- matrix(solve_LSAP(as.matrix(P),  maximum = TRUE))
  P = Matrix::diag(n)
  P = rbind(cbind(Matrix::diag(m), matrix(0, m, n)),
    cbind(matrix(0, n, m), P[corr, ]))
  corr <- cbind(matrix((m+1):totv,  n), matrix(m+corr, n))
  list(corr=corr,  P=P,  D=D)
}
