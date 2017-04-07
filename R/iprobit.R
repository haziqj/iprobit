library(iprior)
library(ggplot2)
library(gridExtra)

ikernL <- function(Xl, kernel = c("Canonical", "FBM,0.5"), interactions = NULL) {
  Hurst <- splitHurst(kernel)  # get the Hurst coefficient
  kernel <- splitKernel(kernel)  # get the kernel
  kernel <- match.arg(kernel, c("Canonical", "FBM"))
  if (kernel == "FBM") kernelFn <- iprior::fnH3
  else kernelFn <- iprior::fnH2
  p <- length(Xl)
  Hl <- NULL

  for (i in 1:p) {
    Hl[[i]] <- kernelFn(as.matrix(Xl[[i]]))
  }
  Hl
}

iprobit <- function(y, ..., kernel = "Canonical", maxit = 100, stop.crit = 1e-3,
                    silent = FALSE, interactions = NULL) {
  # Prepare kernel matrices ----------------------------------------------------
  Xl <- list(...)
  iprobit.kernel <- ikernL(Xl, kernel, NULL)
  H <- iprobit.kernel[[1]]  # CHANGE THIS FOR FUTURE UPDATES. NOW ONLY SUPPORT
                            # SINGLE LAMBDA

  H2 <- H %*% H
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
  n <- length(y)

  # Set up parameter results ---------------------------------------------------
  lower.bound <- lambda <- rep(NA, maxit)
  error.rate <- alpha <- rep(0, maxit)
  w <- matrix(NA, ncol = n, nrow = maxit)
  ystar <- matrix(NA, ncol = n, nrow = maxit)

  # Initialise -----------------------------------------------------------------
  lambda[1] <- abs(rnorm(1))
  lambda2 <- 1
  alpha[1] <- rnorm(1) #0
  w[1, ] <- rnorm(n) #rep(0, n)
  niter <- 1
  lb.const <- (n + 2 - log(n)) / 2 + log(2 * pi)

  for (t in 1:(maxit - 1)) {
    # Update ystar -------------------------------------------------------------
    eta <- as.numeric(alpha[t] + lambda[t] * H %*% w[t, ])
    thing <- rep(NA, n)
    thing1 <- exp(  # phi(eta) / Phi(eta)
      dnorm(eta[y == 1], log = TRUE) - pnorm(eta[y == 1], log.p = TRUE)
    )
    thing0 <- -exp(  # -1 * {phi(eta) / Phi(-eta)}
      dnorm(eta[y == 0], log = TRUE) - pnorm(-eta[y == 0], log.p = TRUE)
    )
    thing[y == 1] <- thing1
    thing[y == 0] <- thing0
    ystar[t + 1, ] <- eta + thing

    # Update w -----------------------------------------------------------------
    A <- lambda2 * H2 + diag(1, n)
    a <- as.numeric(lambda[t] * crossprod(H, ystar[t + 1, ] - alpha[t]))
    eigenA <- eigen(A)
    V <- eigenA$vec
    u <- eigenA$val + 1e-8  # ensure positive eigenvalues
    uinv.Vt <- t(V) / u
    w[t + 1, ] <- as.numeric(crossprod(a, V) %*% uinv.Vt)
    W <- V %*% uinv.Vt + tcrossprod(w[t + 1, ])

    # Update lambda ------------------------------------------------------------
    ct <- sum(H2 * W)
    d <- as.numeric(crossprod(ystar[t + 1, ] - alpha[t], H) %*% w[t + 1, ])
    lambda[t + 1] <- d / ct
    lambda2 <- 1 / ct + (d / ct) ^ 2

    # Update alpha -------------------------------------------------------------
    alpha[t + 1] <- mean(ystar[t + 1, ] - lambda[t + 1] * H %*% w[t + 1, ])

    # Lower bound --------------------------------------------------------------
    lower.bound[t + 1] <- lb.const +
      sum(pnorm(eta[y == 1], log.p = TRUE)) + sum(pnorm(-eta[y == 0], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus + log(ct)) / 2

    # Running fit --------------------------------------------------------------
    tmp <- rep(0, n)
    tmp[ystar[t + 1, ] >= 0] <- 1
    error.rate[t] <- sum(tmp != y) / n * 100

    lb.diff <- abs(lower.bound[t + 1] - lower.bound[t])
    if (!is.na(lb.diff) && (lb.diff < stop.crit)) break
    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, t)
  }
  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar[niter, ], w = w[niter, ], lambda = lambda[niter],
              alpha = alpha[niter], lower.bound = lower.bound, kernel = kernel,
              X = X, y = y, error.rate = error.rate)
  class(res) <- "ipriorProbit"
  res
}

# I-prior probit fitted
fitted.ipriorProbit <- function(x) {
  y.hat <- rep(0, length(x$ystar))
  y.hat[x$ystar >= 0] <- 1
  p.hat <- pnorm(x$ystar)

  list(y = y.hat, prob = p.hat)
}

# I-prior probit predict
predict.ipriorProbit <- function(object, newdata, ...) {
  w <- object$w
  lambda <- object$lambda
  alpha <- object$alpha

  if (object$kernel == "Canonical") H.tilde <- fnH2(object$X, newdata)
  if (object$kernel == "FBM") H.tilde <- fnH3(object$X, newdata)

  ystar.hat <- as.numeric(alpha + lambda * H.tilde %*% w)
  y.hat <- rep(0, nrow(newdata)); y.hat[ystar.hat >= 0] <- 1
  p.hat <- pnorm(ystar.hat)

  list(y = y.hat, prob = p.hat)
}

# I-prior probit plot
plot.ipriorProbit <- function(x, niter.plot = NULL, levels = NULL, ...) {
  if (is.null(niter.plot)) niter.plot <- length(x$lower.bound) - 1
  tmp <- as.factor(x$y)
  if (!is.null(levels)) levels(tmp) <- levels
  lb <- x$lower.bound[2:(niter.plot + 1)]
  error.rate <- x$error.rate[2:(niter.plot + 1)]
  maximin <- max(lb) - min(lb)
  maximin.inv <- 1 / maximin
  error.rate.scaled <- error.rate * maximin + min(lb)
  plot.df1 <- data.frame(Iteration = 1:niter.plot,
                         lower = lb,
                         error = error.rate.scaled)
  plot.df2 <- data.frame(Observation = 1:length(x$ystar),
                         p.hat = fitted(x)$prob,
                         Class = tmp)

  p1 <- ggplot(plot.df1) +
    geom_point(aes(x = Iteration, y = error, col = "Error rate")) +
    geom_line(aes(x = Iteration, y = error, col = "Error rate",
                  linetype = "Error rate")) +
    geom_point(aes(x = Iteration, y = lower, col = "Lower bound")) +
    geom_line(aes(x = Iteration, y = lower, col = "Lower bound",
                  linetype = "Lower bound")) +
    scale_linetype_manual(name = NULL,
                          values = c("Lower bound" = "longdash",
                                     "Error rate" = "solid")) +
    scale_colour_manual(name = NULL,
                        values = c("Lower bound" = "black",
                                   "Error rate" = "lightgoldenrod4")) +
    scale_y_continuous(
      "Lower bound",
      sec.axis = sec_axis(~ (. - min(lb)) * maximin.inv, name = "Error rate")
    ) +
    theme(legend.position = "top")

  p2 <- ggplot(plot.df2, aes(x = Observation, y = p.hat, col = Class)) +
  geom_point() +
  labs(y = "Fitted probabilities")
  #
  # grid.arrange(p1, p2, ncol = 1, nrow = 2, heights = c(6, 4))
  p2
}

print.ipriorProbit <- function(x, newdata = NULL, testdata = NULL) {
  y.hat <- fitted(x)$y
  train.error.rate <- format(round(mean(y.hat != x$y) * 100, 2))
  cat("Training error rate:", train.error.rate, "%")
  if (!is.null(newdata) && !is.null(testdata)) {
    y.hat <- predict(x, newdata)$y
    test.error.rate <- format(round(mean(y.hat != testdata) * 100, 2))
    cat("\nTest error rate:", test.error.rate, "%")
  }
}

# Extract lower bound
logLik.ipriorProbit <- function(x) {
  lb <- x$lower.bound[!is.na(x$lower.bound)]
  cat("Lower bound =", lb[length(lb)])
}