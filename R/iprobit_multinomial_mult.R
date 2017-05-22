data(iris)
# iris <- iris[c(1:3, 51:53, 101:103), ]
y <- iris$Species
y.lev <- levels(y)
X <- iris[, -5]

# MULTIPLE LAMBDA AND ALPHA

n <- length(y)
m <- length(y.lev)
nm <- n * m
p <- ncol(X)
maxit <- 100

# Kernel matrices
H <- iprior::fnH3(X)  # Canonical or FBM
H.sq <- H %*% H

# Initialise
lambda <- rep(1, m)
lambda.sq <- lambda ^ 2
alpha <- rnorm(m)
w <- matrix(0, ncol = m, nrow = n)
f.tmp <- ystar.tmp <- matrix(0, ncol = m, nrow = n)
ystar <- rep(NA, nm)
H.lam <- H.lam.sq <- list(NULL)
for (j in 1:m) {
  H.lam[[j]] <- lambda[j] * H
  H.lam.sq[[j]] <- lambda.sq[j] * H.sq
}
lb <- rep(NA, maxit)
dt <- ct <- logdetA <- rep(NA, m)
W <- list(NULL)
logClb <- rep(NA, n)
lb <- rep(NA, maxit)

for (t in 1:maxit) {
  # Update f
  for (j in 1:m) {
    f.tmp[, j] <- alpha[j] + H.lam[[j]] %*% w[, j]
  }

  # Update ystar
  for (i in 1:n) {
    j <- as.numeric(y[i])
    fi <- f.tmp[i, ]
    fik <- fi[-j]; fij <- fi[j]
    logClb[i] <- logC <- EprodPhiZ(fij - fik, log = TRUE)
    for (k in seq_len(m)[-j]) {
      logD <- log(integrate(
        function(z) {
          fij.minus.fil <- fij - fi[-c(k, j)]
          logPhi.l <- 0
          for (kk in seq_len(length(fij.minus.fil)))
            logPhi.l <- logPhi.l + pnorm(z + fij.minus.fil[kk], log.p = TRUE)
          logphi.k <- dnorm(z + fij - fi[k], log = TRUE)
          exp(logphi.k + logPhi.l) * dnorm(z)
        }, lower = -Inf, upper = Inf
      )$value)
      ystar.tmp[i, k] <- fi[k] - exp(logD - logC)
    }
    ystar.tmp[i, j] <- fi[j] - sum(ystar.tmp[i, -j] - fi[-j])
  }
  ystar <- c(t(ystar.tmp))

  # Update w
  for (j in 1:m) {
    A <- H.lam.sq[[j]] + diag(1, n)
    a <- H.lam[[j]] %*% (ystar.tmp[, j] - alpha[j])
    w[, j] <- solve(A, a)
    W[[j]] <- solve(A) + tcrossprod(w[, j])
    logdetA[j] <- determinant(A)$mod
  }

  # Update lambda
  for (j in 1:m) {
    ct[j] <- sum(H.sq * W[[j]])
    dt[j] <- t(ystar.tmp[, j] - alpha[j]) %*% (H %*% w[, j])
  }
  lambda <- dt / ct
  lambda.sq <- 1 / ct + lambda ^ 2

  # Update H.lam and H.lam.sq
  for (j in 1:m) {
    H.lam[[j]] <- lambda[j] * H
    H.lam.sq[[j]] <- lambda.sq[j] * H.sq
  }

  # Update alpha
  alpha <- apply(ystar.tmp - rep(lambda, each = n) * (H %*% w), 2, mean)

  # Calculate lower bound
  lb.ystar <- sum(logClb)
  lb.w <- 0.5 * (nm - sum(sapply(W, function(x) sum(diag(x)))) - sum(logdetA))
  lb.lambda <- (m / 2) * (1 + log(2 * pi) - mean(log(ct)))
  lb.alpha <- (m / 2) * (1 + log(2 * pi) - log(n))
  lb[t] <- lb.ystar + lb.w + lb.lambda + lb.alpha
}

y.hat <- factor(apply(ystar.tmp, 1, function(x) which(x == max(x))))
levels(y.hat) <- y.lev
probs <- ystar.tmp
for (i in 1:n) {
  for (j in 1:m) {
    probs[i, j] <- exp(
      sum(pnorm((ystar.tmp[i, j] - ystar.tmp[i, (1:m)[-j]]) / sqrt(2),
                log.p = TRUE))
    )
  }
  # probs[i, m] <- 1 - sum(probs[i, 1:(m - 1)])
}
probs <- probs / matrix(rep(apply(probs, 1, sum), m), ncol = m)  # normalise
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
