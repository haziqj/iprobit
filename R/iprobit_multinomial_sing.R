data(iris)
# iris <- iris[c(1:3, 51:53, 101:103), ]
y <- iris$Species
y.lev <- levels(y)
X <- iris[, -5]

# SINGLE LAMBDA AND ALPHA

n <- length(y)
m <- length(y.lev)
nm <- n * m
p <- ncol(X)
maxit <- 100

# Kernel matrices
H <- iprior::fnH2(X)  # Canonical or FBM
H.sq <- H %*% H

# Initialise
lambda <- abs(rnorm(1))
lambda.sq <- lambda ^ 2
alpha <- rnorm(1)
w <- matrix(0, ncol = m, nrow = n)
f.tmp <- ystar.tmp <- matrix(NA, ncol = m, nrow = n)
ystar <- rep(NA, nm)
H.lam <- lambda * H
H.lam.sq <- lambda.sq * H.sq
lb <- rep(NA, maxit)
dt <- ct <- logdetA <- rep(NA, m)
W <- list(NULL)
logClb <- rep(NA, n)
lb <- rep(NA, maxit)

for (t in 1:maxit) {
  # Update f
  f.tmp <- alpha + H.lam %*% w

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
  A <- H.lam.sq + diag(1, n)
  for (j in 1:m) {
    a <- H.lam %*% (ystar.tmp[, j] - alpha)
    w[, j] <- solve(A, a)
    W[[j]] <- solve(A) + tcrossprod(w[, j])
    logdetA[j] <- determinant(A)$mod
  }

  # Update lambda
  for (j in 1:m) {
    ct[j] <- sum(H.sq * W[[j]])
    dt[j] <- t(ystar.tmp[, j] - alpha) %*% (H %*% w[, j])
  }
  lambda <- sum(dt) / sum(ct)
  lambda.sq <- 1 / sum(ct) + lambda ^ 2

  # Update H.lam and H.lam.sq
  H.lam <- lambda * H
  H.lam.sq <- lambda.sq * H.sq

  # Update alpha
  alpha <- mean(ystar.tmp - H.lam %*% w)

  # Calculate lower bound
  lb[t] <- sum(logClb) + nm / 2 - 0.5 * sum(sapply(W, function(x) sum(diag(x)))) -
    0.5 * sum(logdetA) + 0.5 * (1 + log(2 * pi)) - 0.5 * log(sum(ct)) +
    0.5 * (1 + log(2 * pi) - log(nm))
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


