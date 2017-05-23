data(iris)
# iris <- iris[c(1:2, 51:52, 101:102), ]
y <- iris$Species
y.lev <- levels(y)
X <- iris[, -5]

n <- length(y)
m <- length(y.lev)
nm <- n * m
p <- ncol(X)
maxit <- 100

# Expand data ----------------------------------------------------------------
x.expand <- cbind(X, class = y)[rep(1:n, each = m), ]
x.expand <- cbind(x.expand, i = rep(1:n, each = m), j = rep(1:m, times = n))
x.expand <- cbind(x.expand,
                  alpha = as.numeric(x.expand$j == as.numeric(x.expand$class)))
row.names(x.expand) <- NULL

# Kernel matrices ------------------------------------------------------------
H1 <- iprior::fnH2(x.expand[, 1:p])  # main effects (FBM
# H2 <- iprior::fnH1(x.expand$class)  # Pearson
H2 <- iprior::fnH1(as.factor(x.expand$j))
H12 <- H1 * iprior::fnH1(as.factor(x.expand$j))
H2.sq <- H2 %*% H2
H12.sq <- H12 %*% H12
H2H12 <- H2 %*% H12
H12H2 <- t(H2H12)

# Initialise -----------------------------------------------------------------
lambda <- rnorm(2)
lambda.sq <- lambda ^ 2
alpha <- rnorm(1)
# w <- rnorm(nm)
# w <- rep(1, nm)
w <- x.expand$alpha
ystar.tmp <- matrix(NA, ncol = m, nrow = n)
ystar <- rep(NA, nm)
logClb <- rep(NA, n)
H.lam <- lambda[2] * H2 + lambda[1] * H12
H.lam.sq <- lambda.sq[2] * H2.sq + lambda.sq[1] * H12.sq +
  lambda[1] * lambda[2] * (H2H12 + t(H2H12))
lb <- rep(NA, maxit)
niter <- 1

for (t in 1:maxit) {
  # Update ystar -------------------------------------------------------------
  f <- alpha + H.lam %*% w
  for (i in 1:n) {
    j <- as.numeric(y[i])
    fi <- as.numeric(f)[x.expand$i == i]
    fik <- fi[-j]; fij <- fi[j]
    logClb[i] <- logC <- sum(pnorm((fij - fik) / sqrt(2), log.p = TRUE))
    for (k in seq_len(m)[-j]) {
      logD <- log(integrate(
        function(z) {
          logPhi.l <- sum(pnorm(z + fij - fi[-c(k, j)], log.p = TRUE))
          logphi.k <- dnorm(z + fij - fi[k], log = TRUE)
          exp(logphi.k + logPhi.l) * dnorm(z)
        }, lower = -Inf, upper = Inf
      )$value)
      ystar.tmp[i, k] <- fi[k] - exp(logD - logC)
    }
    ystar.tmp[i, j] <- fi[j] - sum(ystar.tmp[i, -j] - fi[-j])
  }
  ystar <- c(t(ystar.tmp))

  # Update w -----------------------------------------------------------------
  A <- H.lam.sq + diag(1, nm)
  a <- H.lam %*% (ystar - alpha)
  logdetA <- determinant(A)$mod
  eigenA <- eigen(A)
  V <- eigenA$vec
  u <- eigenA$val + 1e-8  # ensure positive eigenvalues
  uinv.Vt <- t(V) / u
  w <- as.numeric(crossprod(a, V) %*% uinv.Vt)
  Varw <- V %*% uinv.Vt
  W <- Varw + tcrossprod(w)

  # Update lambda_1 ----------------------------------------------------------
  c1 <- as.numeric(sum(H12.sq * W))
  d1 <- as.numeric(t(ystar - alpha) %*% H12 %*% w - lambda[2] * sum(H2H12 * W))
  lambda[1] <- d1 / c1
  lambda.sq[1] <- 1 / c1 + (d1 / c1) ^ 2

  # Update lambda_2 ----------------------------------------------------------
  c2 <- as.numeric(sum(H2.sq * W))
  d2 <- as.numeric(t(ystar - alpha) %*% H2 %*% w - lambda[1] * sum(H12H2 * W))
  lambda[2] <- d2 / c2
  lambda.sq[2] <- 1 / c2 + (d2 / c2) ^ 2

  # Update H.lam and H.lam.sq ------------------------------------------------
  H.lam <- lambda[2] * H2 + lambda[1] * H12
  H.lam.sq <- lambda.sq[2] * H2.sq + lambda.sq[1] * H12.sq +
    lambda[1] * lambda[2] * (H2H12 + H12H2)

  # Update alpha -------------------------------------------------------------
  alpha <- mean(ystar - H.lam %*% w)

  # Calculate lower bound ----------------------------------------------------
  lb[t + 1] <- 0.5 * (nm - log(nm) + 3 * (1 + log(2 * pi))) -
    0.5 * (logdetA + sum(diag(W)) + log(c1) + log(c2)) + sum(logC)
}


y.hat <- factor(apply(ystar.tmp, 1, function(x) which(x == max(x))))
levels(y.hat) <- y.lev
probs <- ystar.tmp
for (i in 1:n) {
  for (j in 1:(m )) {
    probs[i, j] <- exp(sum(pnorm((ystar.tmp[i, j] -
                                    ystar.tmp[i, (1:m)[-j]]) / sqrt(2), log.p = TRUE)))
  }
  # probs[i, m] <- 1 - sum(probs[i, 1:(m - 1)])
}
# probs <- probs / matrix(rep(apply(probs, 1, sum), m), ncol = m)  # normalise
colnames(probs) <- y.lev
round(probs, 3)

# plot
df.plot <- data.frame(probs, i = 1:n)
df.plot <- reshape2::melt(df.plot, id.vars = "i")
ggplot(df.plot, aes(x = i, y = value)) +
  geom_area(aes(col = variable, fill = variable), position = "stack") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(col = "Class", fill = "Class", x = NULL, y = "Fitted probability") +
  theme_bw()
# ggplot(df.plot, aes(x = i, y = value)) +
#   geom_point(aes(col = variable), position = "stack") +
#   labs(col = "Class", x = NULL, y = "Fitted probability") +
#   coord_cartesian(ylim = c(0, 1)) +
#   facet_grid(variable ~ .) +
#   theme_bw()
