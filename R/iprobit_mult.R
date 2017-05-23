#' @export
iprobit_mult <- function(y, X, kernel = c("Canonical", "FBM"), maxit = 100,
                         stop.crit = 1e-5,  silent = FALSE, alpha0 = NULL,
                         lambda0 = NULL, w0 = NULL, common.intercept = FALSE,
                         common.RKHS.scale = TRUE) {
  n <- length(y)
  y.lev <- levels(y)
  m <- length(y.lev)
  nm <- n * m
  p <- ncol(X)
  kernel <- match.arg(kernel, c("Canonical", "FBM"))
  if (kernel == "FBM") H <- iprior::fnH3(X)
  else H <- iprior::fnH2(X)
  H.sq <- H %*% H

  # Initialise
  alpha <- rnorm(m)
  lambda <- abs(rnorm(m))
  if (isTRUE(common.intercept)) alpha <- rep(rnorm(1), m)
  if (isTRUE(common.RKHS.scale)) lambda <- rep(rnorm(1), m)
  lambda.sq <- lambda ^ 2
  H.lam <- H.lam.sq <- list(NULL)
  for (j in 1:m) {
    H.lam[[j]] <- lambda[j] * H
    H.lam.sq[[j]] <- lambda.sq[j] * H.sq
  }
  w <- f.tmp <- ystar.tmp <- matrix(0, ncol = m, nrow = n)
  dt <- ct <- logdetA <- rep(NA, m)
  lb <- rep(NA, maxit)
  W <- list(NULL)
  logClb <- rep(NA, n)
  niter <- 1

  # Supplied starting values
  if (!is.null(alpha0)) alpha <- alpha0
  if (!is.null(lambda0)) lambda <- lambda0
  if (!is.null(w0)) w <- w0

  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
  start.time <- Sys.time()
  for (t in 1:(maxit - 1)) {
    # Update f
    f.tmp <- rep(alpha, each = n) + rep(lambda, each = n) * (H %*% w)

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
    if (isTRUE(common.RKHS.scale)) {
      lambda <- rep(sum(dt) / sum(ct), m)
      lambda.sq <- rep(1 / sum(ct) + lambda ^ 2, m)
    } else {
      lambda <- dt / ct
      lambda.sq <- 1 / ct + lambda ^ 2
    }

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
    if (isTRUE(common.RKHS.scale))
      lb.lambda <- 0.5 * (1 + log(2 * pi) - log(sum(ct)))
    else
      lb.lambda <- (m / 2) * (1 + log(2 * pi) - mean(log(ct)))
    if (isTRUE(common.RKHS.scale))
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
  time.taken <- end.time - start.time

  # Clean up and close
  if (isTRUE(common.RKHS.scale)) lambda <- lambda[1]
  if (isTRUE(common.intercept)) alpha <- alpha[1]
  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar.tmp, w = w, lambda = lambda, alpha = alpha,
              lower.bound = lb, time = time.taken, niter = niter, y = y, X = X,
              kernel = kernel, maxit = maxit)
  class(res) <- "iprobitMult"
  res
}

#' @export
print.iprobitMult <- function(x) {
  Intercept <- lambda <- 0
  m <- ncol(x$w)
  Intercept[1:m] <- x$alpha
  lambda[1:m] <- x$lambda
  theta <- rbind(Intercept, lambda)
  colnames(theta) <- paste0("Class = ", seq_len(m))
  # y.hat <- fitted(x)$y

  cat("Lower bound value = ", x$lower.bound[x$niter], "\n")
  # cat("Training error rate = ", mean(y.hat != x$y), "%\n")
  cat("Iterations = ", x$niter, "\n\n")
  print(round(theta, 5))
}







