sgm <- function(A, B, m, start, iteration) {
  ## Seeds are assumed to be vertices 1:m in both graphs
  totv <- ncol(A)
  n <- totv - m
  if (m != 0) {
    A12 <- rbind(A[1:m, (m+1):(m+n)])
    A21 <- cbind(A[(m+1):(m+n), 1:m])
    B12 <- rbind(B[1:m, (m+1):(m+n)])
    B21 <- cbind(B[(m+1):(m+n), 1:m])
  }
  if ( m==0 ) {
    A12 <- matrix(0, n, n)
    A21 <- matrix(0, n, n)
    B12 <- matrix(0, n, n)
    B21 <- matrix(0, n, n)
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
    x <-  A21 %*% t(B21)
    y <-  t(A12) %*% B12
    z <-  A22 %*% P %*% t(B22)
    w <-  t(A22) %*% P %*% B22
    Grad <- x + y + z + w
    ind <- matrix(solve_LSAP(Grad, maximum = TRUE))
    T <- diag(n)
    T <- T[ind, ]
    wt <- t(A22) %*% T %*% B22
    c <- sum(diag(w %*% t(P)))
    d <- sum(diag(wt %*% t(P)))+sum(diag(w %*% t(T)))
    e <- sum(diag(wt %*% t(T)))
    u <- sum(diag(t(P) %*% x+t(P) %*% y))
    v <- sum(diag(t(T) %*% x+t(T) %*% y))
    if( c-d+e == 0 && d-2*e+u-v == 0) {
      alpha <- 0
    } else {
      alpha <- -(d-2*e+u-v) / (2*(c-d+e))}
    f0 <- 0
    f1 <-  c-e+u-v
    falpha <- (c-d+e) * alpha^2 + (d-2*e+u-v) * alpha
    if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1) {
      P <-  alpha*P + (1-alpha) * T
    } else if(f0 > f1) {
      P <- T
    } else {
      toggle <- 0}
  }
  D <- P
  corr <- matrix(solve_LSAP(P,  maximum = TRUE))
  P = diag(n)
  P = rbind(cbind(diag(m), matrix(0, m, n)), cbind(matrix(0, n, m), P[corr, ]))
  corr <- cbind(matrix((m+1):totv,  n), matrix(m+corr, n))
  return(list(corr=corr,  P=P,  D=D))
}
