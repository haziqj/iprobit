#' @export
iprobit_mult <- function(ipriorKernel, maxit = 100, stop.crit = 1e-5,
                         silent = FALSE, alpha0 = NULL, lambda0 = NULL,
                         w0 = NULL, common.intercept = FALSE,
                         common.RKHS.scale = FALSE) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(ipriorKernel, iprobit.env)
  list2env(BlockBstuff, iprobit.env)
  list2env(model, iprobit.env)
  environment(BlockB) <- iprobit.env
  environment(lambdaExpand_mult) <- iprobit.env
  environment(HlamFn_mult) <- iprobit.env
  environment(HlamsqFn_mult) <- iprobit.env
  y <- Y
  m <- length(y.levels)
  nm <- n * m

  # Initialise -----------------------------------------------------------------
  if (is.null(alpha0)) {
    if (isTRUE(common.intercept)) alpha0 <- rep(rnorm(1), m)
    else alpha0 <- rnorm(m)
  }
  if (is.null(lambda0)) {
    if (isTRUE(common.RKHS.scale)) lambda0 <- rep(abs(rnorm(l)), m)
    else lambda0 <- abs(rnorm(l * m))
  }
  if (is.null(w0)) w0 <- matrix(0, ncol = m, nrow = n)
  alpha <- alpha0
  lambda <- ct <- dt <- matrix(lambda0, ncol = m, nrow = l)
  lambda.sq <- lambda ^ 2
  lambdaExpand_mult(env = iprobit.env)
  HlamFn_mult(env = iprobit.env)
  HlamsqFn_mult(env = iprobit.env)
  w <- f.tmp <- ystar <- w0
  logdetA <- rep(NA, m)
  lb <- rep(NA, maxit)
  W <- list(NULL)
  logClb <- rep(NA, n)
  niter <- 1

  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
  start.time <- Sys.time()
  for (t in 1:(maxit - 1)) {
    # Update f
    f.tmp <- rep(alpha, each = n) + mapply("%*%", Hlam.mat, split(w, col(w)))

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
        ystar[i, k] <- fi[k] - exp(logD - logC)
      }
      ystar[i, j] <- fi[j] - sum(ystar[i, -j] - fi[-j])
    }

    # Update w
    for (j in 1:m) {
      A <- Hlam.matsq[[j]] + diag(1, n)
      a <- Hlam.mat[[j]] %*% (ystar[, j] - alpha[j])
      w[, j] <- solve(A, a)
      W[[j]] <- solve(A) + tcrossprod(w[, j])
      logdetA[j] <- determinant(A)$mod
    }

    # Update lambda
    for (k in 1:l) {
      for (j in 1:m) {
        lambdaExpand_mult(env = iprobit.env)
        BlockB(k, lambda[, j])
        ct[k, j] <- sum(Psql[[k]] * W[[j]])
        dt[k, j] <- as.numeric(
          crossprod(ystar[, j] - alpha[j], Pl[[k]]) %*% w[, j] -
            sum(Sl[[k]] * W[[j]]) / 2
        )
      }
      if (isTRUE(common.RKHS.scale)) {
        lambda[k, ] <- rep(sum(dt[k, ]) / sum(ct[k, ]), m)
        lambda.sq[k, ] <- rep(1 / sum(ct[k, ]) + lambda[k, 1] ^ 2, m)
      } else {
        lambda[k, ] <- dt[k, ] / ct[k, ]
        lambda.sq[k, ] <- 1 / ct[k, ] + lambda[k, ] ^ 2
      }
    }

    # Update H.lam and H.lam.sq
    lambdaExpand_mult(env = iprobit.env)
    HlamFn_mult(env = iprobit.env)
    HlamsqFn_mult(env = iprobit.env)

    # Update alpha
    alpha <- apply(ystar - mapply("%*%", Hlam.mat, split(w, col(w))), 2, mean)
    if (isTRUE(common.intercept)) alpha <- rep(mean(alpha), m)

    # Calculate lower bound
    lb.ystar <- sum(logClb)
    lb.w <- 0.5 * (nm - sum(sapply(W, function(x) sum(diag(x)))) - sum(logdetA))
    if (isTRUE(common.RKHS.scale))
      lb.lambda <- (l / 2) * (1 + log(2 * pi)) - sum(log(apply(ct, 1, sum))) / 2
    else
      lb.lambda <- (m / 2) * (l * (1 + log(2 * pi)) - sum(log(ct)) / m)
    if (isTRUE(common.intercept))
      lb.alpha <- 0.5 * (1 + log(2 * pi) - log(nm))
    else
      lb.alpha <- (m / 2) * (1 + log(2 * pi) - log(n))
    lb[t + 1] <- lb.ystar + lb.w + lb.lambda + lb.alpha

    # Lower bound difference
    lb.diff <- abs(lb[t + 1] - lb[t])
    if (!is.na(lb.diff) && (lb.diff < stop.crit)) break
    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, t)
  }
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Calculate standard errors from posterior variance --------------------------
  if (isTRUE(common.intercept)) se.alpha <- sqrt(1 / nm)
  else se.alpha <- sqrt(1 / n)
  if (isTRUE(common.RKHS.scale))
    se.lambda <- matrix(sqrt(1 / apply(ct, 1, sum)), ncol = m, nrow = l)
  else
    se.lambda <- matrix(sqrt(1 / ct[1:l, ]), ncol = m, nrow = l)
  se.ystar <- NA #iprobitSE(y = y, eta = eta, thing1 = thing1, thing0 = thing0)

  # Clean up and close
  lambda <- matrix(lambda[1:l, ], ncol = m, nrow = l)
  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar, w = w, lambda = lambda, alpha = alpha,
              lower.bound = lb, ipriorKernel = ipriorKernel, se.alpha = se.alpha,
              se.lambda = se.lambda, se.ystar = se.ystar,
              y.levels = y.levels, start.time = start.time,
              end.time = end.time, time = time.taken, call = match.call(),
              stop.crit = stop.crit, niter = niter, maxit = maxit)
  class(res) <- c("iprobitMod", "iprobitMod_mult")
  res
}

