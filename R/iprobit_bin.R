#' @export
iprobit_bin <- function(ipriorKernel, maxit = 200, stop.crit = 1e-5,
                        silent = FALSE, alpha0 = NULL, lambda0 = NULL,
                        w0 = NULL) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(ipriorKernel, iprobit.env)
  list2env(BlockBstuff, iprobit.env)
  list2env(model, iprobit.env)
  environment(BlockB) <- iprobit.env
  environment(.lambdaExpand) <- iprobit.env
  environment(HlamFn) <- iprobit.env
  environment(HlamsqFn) <- iprobit.env
  y <- Y

  # Initialise -----------------------------------------------------------------
  if (is.null(lambda0)) lambda0 <- abs(rnorm(l))
  if (is.null(alpha0)) alpha0 <- rnorm(1)
  if (is.null(w0)) w0 <- rep(0, n)
  lambda <- ct <- lambda0
  lambda.sq <- lambda0 ^ 2
  .lambdaExpand(env = iprobit.env)
  HlamFn(iprobit.env)
  HlamsqFn(iprobit.env)
  alpha <- alpha0
  w <- w0
  niter <- 1
  lower.bound <- rep(NA, maxit)
  lb.const <- (n + 1 + l - log(n) + (l + 1) * log(2 * pi)) / 2

  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
  start.time <- Sys.time()
  for (t in 1:(maxit - 1)) {
    # Update ystar -------------------------------------------------------------
    eta <- as.numeric(alpha + Hlam.mat %*% w)
    thing <- rep(NA, n)
    thing1 <- exp(  # phi(eta) / Phi(eta)
      dnorm(eta[y == 2], log = TRUE) - pnorm(eta[y == 2], log.p = TRUE)
    )
    thing0 <- -exp(  # -1 * {phi(eta) / Phi(-eta)}
      dnorm(eta[y == 1], log = TRUE) - pnorm(-eta[y == 1], log.p = TRUE)
    )
    thing[y == 2] <- thing1
    thing[y == 1] <- thing0
    ystar <- eta + thing

    # Update w -----------------------------------------------------------------
    A <- Hlam.matsq + diag(1, n)
    a <- as.numeric(crossprod(Hlam.mat, ystar - alpha))
    eigenA <- iprior::eigenCpp(A)
    V <- eigenA$vec
    u <- eigenA$val + 1e-8  # ensure positive eigenvalues
    uinv.Vt <- t(V) / u
    w <- as.numeric(crossprod(a, V) %*% uinv.Vt)
    Varw <- iprior::fastVDiag(V, 1 / u)  # V %*% uinv.Vt
    W <- Varw + tcrossprod(w)

    # Update lambda ------------------------------------------------------------
    for (k in 1:l) {
      .lambdaExpand(env = iprobit.env)
      BlockB(k)
      ct[k] <- sum(Psql[[k]] * W)
      d <- as.numeric(
        crossprod(ystar - alpha, Pl[[k]]) %*% w - sum(Sl[[k]] * W) / 2
      )
      lambda[k] <- d / ct[k]
      lambda.sq[k] <- 1 / ct[k] + (d / ct[k]) ^ 2
    }
    .lambdaExpand(env = iprobit.env)
    HlamFn(iprobit.env)
    HlamsqFn(iprobit.env)

    # Update alpha -------------------------------------------------------------
    alpha <- mean(ystar - Hlam.mat %*% w)

    # Lower bound --------------------------------------------------------------
    lower.bound[t + 1] <- lb.const +
      sum(pnorm(eta[y == 2], log.p = TRUE)) +
      sum(pnorm(-eta[y == 1], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus + sum(log(ct))) / 2

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
  se.ystar <- NA #iprobitSE(y = y, eta = eta, thing1 = thing1, thing0 = thing0)

  # Close function -------------------------------------------------------------
  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar, w = w, lambda = lambda[1:l], alpha = alpha,
              lower.bound = lower.bound, kernel = kernel, ipriorKernel = ipriorKernel,
              se = c(se.alpha, se.lambda), se.ystar = se.ystar,
              y.levels = y.levels, Varw = Varw, start.time = start.time,
              end.time = end.time, time = time.taken, call = match.call(),
              stop.crit = stop.crit, niter = niter, maxit = maxit)
  class(res) <- c("iprobitMod", "iprobitMod_bin")
  res
}
