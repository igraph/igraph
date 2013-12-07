
adjcorr <- function(A, P, corr, permutation) {
  ## input is A modelled from a random binomial graph with A_{i,j}
  ## distributed Bin(P_{i,j}) for example, if P=.5*matrix(1,n,n) then
  ## A is ER(n,0.5) output is B which is adjacency matrix with
  ## correlation corr element-wise to A the labels of B are then
  ## permuted via permutation
  Q <- P + corr * (1 - P)
  n <- nrow(A)
  B <- A
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (A[i, j] == 1 && runif(1) > Q[i, j]) {
        B[i, j] <- 0
        B[j, i] <- 0
      } else if (A[i, j] == 0 && runif(1) < ((1 - Q[i, j]) * (P[i, j]/(1 - 
                    P[i, j])))) {
        B[i, j] <- 1
        B[j, i] <- 1
      }
    }
  }
  P <- diag(n)
  P <- P[permutation, ]
  B <- P %*% B %*% t(P)
  B
}

sgm <- function(A, B, m, start, iteration) {
    # seeds are assumed to be vertices 1:m in both graphs
    totv <- ncol(A)
    n <- totv - m
    A12 <- A[1:m, (m + 1):(m + n)]
    A21 <- as.matrix(A[(m + 1):(m + n), 1:m])
    A22 <- A[(m + 1):(m + n), (m + 1):(m + n)]
    B12 <- B[1:m, (m + 1):(m + n)]
    B21 <- as.matrix(B[(m + 1):(m + n), 1:m])
    B22 <- B[(m + 1):(m + n), (m + 1):(m + n)]

    patience <- iteration
    tol <- 0.99
    P <- start
    toggle <- 1
    iter <- 0
    while (toggle == 1 & iter < patience) {
        iter <- iter + 1
        Grad <- 2 * A22 %*% P %*% t(B22) + 2 * A21 %*% t(B21)
        ind <- matrix(clue::solve_LSAP(Grad, maximum = TRUE))
        T <- diag(n)
        T <- T[ind, ]
        c <- sum(diag(t(A22) %*% P %*% B22 %*% t(P)))
        d <- sum(diag(t(A22) %*% T %*% B22 %*% t(P))) + sum(diag(t(A22) %*% 
            P %*% B22 %*% t(T)))
        e <- sum(diag(t(A22) %*% T %*% B22 %*% t(T)))
        u <- 2 * sum(diag(t(P) %*% A21 %*% t(B21)))
        v <- 2 * sum(diag(t(T) %*% A21 %*% t(B21)))
        if (c - d + e == 0 && d - 2 * e + u - v == 0) {
            alpha <- 0
        } else {
            alpha <- -(d - 2 * e + u - v)/(2 * (c - d + e))
        }
        f0 <- 0
        f1 <- c - e + u - v
        falpha <- (c - d + e) * alpha^2 + (d - 2 * e + u - v) * alpha
        if (alpha < tol && alpha > 0 && falpha > f0 && falpha > f1) {
            P <- alpha * P + (1 - alpha) * T
        } else if (f0 > f1) {
            P <- T
        } else {
            toggle <- 0
        }
    }
    corr <- matrix(clue::solve_LSAP(P, maximum = TRUE))
    corr <- cbind(matrix((m + 1):totv, n), matrix(m + corr, n))
    return(corr)
}
