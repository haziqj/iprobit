#' @export
iprobit_bin <- function(y, ..., kernel = "Canonical", maxit = 100,
                        stop.crit = 1e-5, silent = FALSE, interactions = NULL,
                        alpha0 = rnorm(1), lambda0 = abs(rnorm(1)),
                        w0 = rep(0, n)) {
  y.tmp <- checkLevels(y)
  y <- y.tmp$y
  y.levels <- y.tmp$levels

  # Prepare kernel matrices ----------------------------------------------------
  Xl <- list(...)
  iprobit.kernel <- ikernL(Xl, NULL, kernel, NULL)
  H <- iprobit.kernel[[1]]  # CHANGE THIS FOR FUTURE UPDATES. NOW ONLY SUPPORT
                            # SINGLE LAMBDA.
  H.sq <- H %*% H

  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
  n <- length(y)

  # Initialise -----------------------------------------------------------------
  lambda <- lambda0
  lambda.sq <- lambda ^ 2
  alpha <- alpha0
  w <- w0
  niter <- 1
  lower.bound <- rep(NA, maxit)
  lb.const <- (n + 2 - log(n)) / 2 + log(2 * pi)

  start.time <- Sys.time()
  for (t in 1:(maxit - 1)) {
    # Update ystar -------------------------------------------------------------
    eta <- as.numeric(alpha + lambda * H %*% w)
    thing <- rep(NA, n)
    thing1 <- exp(  # phi(eta) / Phi(eta)
      dnorm(eta[y == 1], log = TRUE) - pnorm(eta[y == 1], log.p = TRUE)
    )
    thing0 <- -exp(  # -1 * {phi(eta) / Phi(-eta)}
      dnorm(eta[y == 0], log = TRUE) - pnorm(-eta[y == 0], log.p = TRUE)
    )
    thing[y == 1] <- thing1
    thing[y == 0] <- thing0
    ystar <- eta + thing

    # Update w -----------------------------------------------------------------
    A <- lambda.sq * H.sq + diag(1, n)
    a <- as.numeric(lambda * crossprod(H, ystar - alpha))
    eigenA <- eigen(A)
    V <- eigenA$vec
    u <- eigenA$val + 1e-8  # ensure positive eigenvalues
    uinv.Vt <- t(V) / u
    w <- as.numeric(crossprod(a, V) %*% uinv.Vt)
    Varw <- V %*% uinv.Vt
    W <- Varw + tcrossprod(w)

    # Update lambda ------------------------------------------------------------
    ct <- sum(H.sq * W)
    d <- as.numeric(crossprod(ystar - alpha, H) %*% w)
    lambda <- d / ct
    lambda.sq <- 1 / ct + (d / ct) ^ 2

    # Update alpha -------------------------------------------------------------
    alpha <- mean(ystar - lambda * H %*% w)

    # Lower bound --------------------------------------------------------------
    lower.bound[t + 1] <- lb.const +
      sum(pnorm(eta[y == 1], log.p = TRUE)) +
      sum(pnorm(-eta[y == 0], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus + log(ct)) / 2

    lb.diff <- abs(lower.bound[t + 1] - lower.bound[t])
    if (!is.na(lb.diff) && (lb.diff < stop.crit)) break
    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, t)
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  # Calculate standard errors from posterior variance --------------------------
  se.lambda <- sqrt(1 / ct)
  se.alpha <- sqrt(1 / n)
  se.ystar <- iprobitSE(y = y, eta = eta, thing1 = thing1, thing0 = thing0)

  # Close function -------------------------------------------------------------
  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar, w = w, lambda = lambda, alpha = alpha,
              lower.bound = lower.bound, kernel = kernel, X = Xl, y = y,
              se = c(se.alpha, se.lambda), se.ystar = se.ystar,
              y.levels = y.levels, Varw = Varw, start.time = start.time,
              end.time = end.time, time = time.taken, call = match.call(),
              stop.crit = stop.crit, niter = niter, maxit = maxit)
  class(res) <- c("iprobitMod", "iprobitMod_bin")
  res
}
