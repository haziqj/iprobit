#' @export
iprobit2 <- function(y, ..., kernel = "Canonical", maxit = 10000, stop.crit = 1e-5,
                     silent = FALSE, interactions = TRUE, alpha0 = rnorm(1),
                     lambda0 = abs(rnorm(2)), w0 = rep(0, n)) {
  # A temporary function to do iprobit models with two covariates and
  # interactions.

  y.tmp <- checkLevels(y)
  y <- y.tmp$y
  y.levels <- y.tmp$levels

  # Prepare kernel matrices ----------------------------------------------------
  Xl <- list(...)
  iprobit.kernel <- ikernL(Xl, NULL, kernel, NULL)
  H1 <- iprobit.kernel[[1]]
  H2 <- iprobit.kernel[[2]]
  H1sq <- H1 %*% H1
  H2sq <- H2 %*% H2
  H1.H2 <- H1 %*% H2
  if (isTRUE(interactions)) {
    H12 <- H1 * H2
    H12sq <- H12 %*% H12
    H1.H12 <- H1 %*% H12
    H2.H12 <- H2 %*% H12
  }

  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
  n <- length(y)

  # Set up parameter results ---------------------------------------------------
  lower.bound <- rep(NA, maxit)
  error.rate <- alpha <- rep(0, maxit)
  lambda <- matrix(NA, ncol = 2, nrow = maxit)
  w <- matrix(NA, ncol = n, nrow = maxit)
  ystar <- matrix(NA, ncol = n, nrow = maxit)

  # Initialise -----------------------------------------------------------------
  lambda[1, ] <- lambda0
  lambdasq <- lambda[1, ] ^ 2
  Hlam.mat <- lambda[1, 1] * H1 + lambda[1, 2] * H2
  if (isTRUE(interactions))
    Hlam.mat <- Hlam.mat + lambda[1, 1] * lambda[1, 2] * H12
  alpha[1] <- alpha0
  w[1, ] <- w0
  niter <- 1
  lb.const <- (n + 2 - log(n)) / 2 + log(2 * pi) + (1 / 2) * (1 + log(2 * pi))

  start.time <- Sys.time()
  for (t in 1:(maxit - 1)) {
    # Update ystar -------------------------------------------------------------
    eta <- as.numeric(alpha[t] + Hlam.mat %*% w[t, ])
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
    Hlamsq.mat <- lambdasq[1] * H1sq + lambdasq[2] * H2sq +
      lambda[t, 1] * lambda[t, 2] * (H1.H2 + t(H1.H2))
    if (isTRUE(interactions))
      Hlamsq.mat <- Hlamsq.mat + lambdasq[1] * lambdasq[2] * H12sq +
      lambdasq[1] * lambda[t, 2] * (H1.H12 + t(H1.H12)) +
      lambda[t, 1] * lambdasq[2] * (H2.H12 + t(H2.H12))
    A <- Hlamsq.mat + diag(1, n)
    a <- as.numeric(crossprod(Hlam.mat, ystar[t + 1, ] - alpha[t]))
    eigenA <- eigen(A)
    V <- eigenA$vec
    u <- eigenA$val + 1e-8  # ensure positive eigenvalues
    uinv.Vt <- t(V) / u
    w[t + 1, ] <- as.numeric(crossprod(a, V) %*% uinv.Vt)
    Varw <- V %*% uinv.Vt
    W <- Varw + tcrossprod(w[t + 1, ])

    # Update lambda[1] ---------------------------------------------------------
    P1 <- H1
    P1.H2 <- H1.H2
    S1 <- 0
    ct1 <- sum(H1sq * W)
    if (isTRUE(interactions)) {
      P1 <- P1 + lambda[t, 2] * H12
      P1.H2 <- P1.H2 + lambda[t, 2] * t(H2.H12)
      S1 <- lambda[t, 2] * (P1.H2 + t(P1.H2))
      ct1 <- ct1 + lambda[t, 2] ^ 2 * sum(H12sq * W) +
        lambda[t, 2] * (sum(H1.H12 * W) + sum(t(H1.H12) * W))
    }
    d1 <- crossprod(ystar[t + 1, ] - alpha[t], P1) %*% w[t + 1, ] - sum(S1 * W) / 2
    lambda[t + 1, 1] <- d1 / ct1
    lambdasq[1] <- 1 / ct1 + (d1 / ct1) ^ 2

    # Update lambda[2] ---------------------------------------------------------
    P2 <- H2
    P2.H1 <- t(H1.H2)
    S2 <- 0
    ct2 <- sum(H2sq * W)
    if (isTRUE(interactions)) {
      P2 <- P2 + lambda[t + 1, 1] * H12
      P2.H1 <- P2.H1 + lambda[t + 1, 1] * t(H1.H12)
      S2 <- lambda[t + 1, 1] * (P2.H1 + t(P2.H1))
      ct2 <- ct2 + lambda[t + 1, 1] ^ 2 * sum(H12sq * W) +
        lambda[t + 1, 1] * (sum(H2.H12 * W) + sum(t(H2.H12) * W))
    }
    d2 <- crossprod(ystar[t + 1, ] - alpha[t], P2) %*% w[t + 1, ] - sum(S2 * W) / 2
    lambda[t + 1, 2] <- d2 / ct2
    lambdasq[2] <- 1 / ct2 + (d2 / ct2) ^ 2

    # Update alpha -------------------------------------------------------------
    Hlam.mat <- lambda[t + 1, 1] * H1 + lambda[t + 1, 2] * H2
    if (isTRUE(interactions))
      Hlam.mat <- Hlam.mat + lambda[t + 1, 1] * lambda[t + 1, 2] * H12
    alpha[t + 1] <- mean(ystar[t + 1, ] - Hlam.mat %*% w[t + 1, ])

    # Lower bound --------------------------------------------------------------
    lower.bound[t + 1] <- lb.const +
      sum(pnorm(eta[y == 1], log.p = TRUE)) + sum(pnorm(-eta[y == 0], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus + log(ct1) + log(ct2)) / 2

    # Running fit --------------------------------------------------------------
    tmp <- rep(0, n)
    tmp[ystar[t + 1, ] >= 0] <- 1
    error.rate[t] <- sum(tmp != y) / n * 100

    lb.diff <- abs(lower.bound[t + 1] - lower.bound[t])
    if (!is.na(lb.diff) && (lb.diff < stop.crit)) break
    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, t)
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  # Calculate standard errors from posterior variance --------------------------
  se.lambda <- c(sqrt(1 / ct1), sqrt(1 / ct2))
  se.alpha <- sqrt(1 / n)
  se.ystar <- iprobitSE(y = y, eta = eta, thing1 = thing1, thing0 = thing0)

  # Close function -------------------------------------------------------------

  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar[niter, ], w = w[niter, ], lambda = lambda[niter, ],
              alpha = alpha[niter], lower.bound = lower.bound, kernel = kernel,
              X = Xl[[1]], y = y, error.rate = error.rate,
              se = c(se.alpha, se.lambda),
              se.ystar = se.ystar, y.levels = y.levels, Varw = Varw,
              start.time = start.time, end.time = end.time, time = time.taken,
              call = match.call(), stop.crit = stop.crit, niter = niter,
              maxit = maxit)
  class(res) <- "ipriorProbit"
  res
}
