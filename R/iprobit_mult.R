data(iris)
# iris <- iris[c(1:3, 51:53, 101:103), ]
y <- iris$Species
y.lev <- levels(y)
X <- iris[, -5]

n <- length(y)
m <- length(y.lev)
nm <- n * m
p <- ncol(X)
maxit <- 10

# Expand data
x.expand <- cbind(X,class = y)[rep(1:n, each = m), ]
x.expand <- cbind(x.expand, i = rep(1:n, each = m), j = rep(1:m, times = n))
x.expand <- cbind(x.expand,
                  alpha = as.numeric(x.expand$j == as.numeric(x.expand$class)))
row.names(x.expand) <- NULL

# Kernel matrices
H1 <- iprior::fnH2(x.expand[, -(p + 1)])
H2 <- iprior::fnH1(x.expand[, p + 1])
H12 <- H1 * H2
H2.sq <- H2 %*% H2
H12.sq <- H12 %*% H12
H2H12 <- H2 %*% H12
trH2sq <- sum(diag(H2.sq))
trH12sq <- sum(diag(H12.sq))

# Initialise
lambda <- abs(rnorm(2))
lambda.sq <- lambda ^ 2
alpha <- rnorm(1)  #x.expand$alpha
w <- rep(0, nm)
ystar.tmp <- matrix(NA, ncol = m, nrow = n)
ystar <- rep(NA, nm)
H.lam <- lambda[2] * H2 + lambda[1] * lambda[2] * H12
H.lam.sq <- lambda.sq[2] * H2.sq + lambda.sq[1] * lambda.sq[2] * H12.sq +
  lambda[1] * lambda.sq[2] * (H2H12 + t(H2H12))
lb <- rep(NA, maxit)

for (t in 1:maxit) {
  # Update ystar
  f <- alpha + H.lam %*% w
  for (i in 1:n) {
    j <- as.numeric(y[i])
    fi <- as.numeric(f)[x.expand$i == i]
    fik <- fi[-j]; fij <- fi[j]
    logC <- sum(pnorm((fij - fik) / sqrt(2), log.p = TRUE))
    for (k in seq_len(m)[-j]) {
      logD <- log(integrate(
        f = function(z) {
          logPhi.l <- sum(pnorm(z + fij - fi[-c(k, j)], log.p = TRUE))
          exp(dnorm(z + fij - fi[k], log = TRUE) + logPhi.l + dnorm(z, log = TRUE))
        }, lower = -Inf, upper = Inf
      )$value)
      ystar.tmp[i, k] <- fi[k] - exp(logD - logC)
    }
    ystar.tmp[i, j] <- fi[j] - sum(ystar.tmp[i, -j] - fi[-j])
  }
  ystar <- c(t(ystar.tmp))

  # Update w
  A <- H.lam.sq + diag(1, nm)
  a <- H.lam %*% (ystar - alpha)
  w <- solve(A, a)
  W <- solve(A) + tcrossprod(w)
  logdetA <- determinant(A)$mod

  # Update lambda_1
  c1 <- as.numeric(lambda[2] * sum(H12.sq * W))
  d1 <- as.numeric(
    lambda[2] * t(ystar - alpha) %*% H12 %*% w - lambda.sq[2] * sum(H2H12 * W)
  )
  lambda[1] <- d1 / c1
  lambda.sq[1] <- 1 / c1 + (d1 / c1) ^ 2

  # Update lambda_2
  c2 <- as.numeric(
    sum(H2.sq * W) + lambda.sq[1] * sum(H12.sq * W) +
      lambda[1] * (sum(H2H12 * W) + sum(t(H2H12 * W)))
  )
  d2 <- as.numeric(t(ystar - alpha) %*% (H2 + lambda[1] * H12) %*% w)
  lambda[2] <- d2 / c2
  lambda.sq[2] <- 1 / c2 + (d2 / c2) ^ 2

  # Update H.lam and H.lam.sq
  H.lam <- lambda[2] * H2 + lambda[1] * lambda[2] * H12
  H.lam.sq <- lambda.sq[2] * H2.sq + lambda.sq[1] * lambda.sq[2] * H12.sq +
    lambda[1] * lambda.sq[2] * (H2H12 + t(H2H12))

  # Update alpha
  alpha <- mean(ystar - H.lam %*% w)

  # Calculate lower bound
  lb[t] <- 0.5 * (nm - log(nm) + 3 * (1 + log(2 * pi))) -
    0.5 * (logdetA + sum(diag(W)) + log(c1) + log(c2)) + sum(logC)
}


