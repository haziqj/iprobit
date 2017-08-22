################################################################################
#
#   iprobit: Binary and Multinomial Probit Regression with I-priors
#   Copyright (C) 2017  Haziq Jamil
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

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
  environment(loop_logical) <- iprobit.env
  y <- Y
  m <- length(y.levels)
  nm <- n * m
  maxit <- max(1, maxit)

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
  alpha <- rep(1,m); alpha[1:m] <- alpha0
  lambda <- ct <- dt <- matrix(lambda0, ncol = m, nrow = l)
  lambda.sq <- lambda ^ 2
  lambdaExpand_mult(env = iprobit.env)
  HlamFn_mult(env = iprobit.env)
  HlamsqFn_mult(env = iprobit.env)
  w <- f.tmp <- ystar <- w0
  W <- list(NULL)
  logdetA <- rep(NA, m)
  w.ind <- seq_len(
    ifelse(isTRUE(common.intercept) && isTRUE(common.RKHS.scale), 1, m)
  )
  niter <- 0
  lb <- error.rates <- brier.scores <- rep(NA, maxit)
  logClb <- rep(NA, n)

  # The variational EM loop ----------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (loop_logical()) {  # see loop_logical() function in iprobit_helper.R
    # Update f -----------------------------------------------------------------
    f.tmp <- rep(alpha, each = n) + mapply("%*%", Hlam.mat, split(w, col(w)))

    # Update ystar -------------------------------------------------------------
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

    # Update w -----------------------------------------------------------------
    for (j in w.ind) {
      A <- Hlam.matsq[[j]] + diag(1, n)
      a <- as.numeric(crossprod(Hlam.mat[[j]], ystar[, j] - alpha[j]))
      if (!isNystrom(ipriorKernel)) {
        eigenA <- iprior::eigenCpp(A)
        V <- eigenA$vec
        u <- eigenA$val + 1e-8  # ensure positive eigenvalues
        uinv.Vt <- t(V) / u
        w[, j] <- as.numeric(crossprod(a, V) %*% uinv.Vt)
        Varw <- iprior::fastVDiag(V, 1 / u)  # V %*% uinv.Vt
        W[[j]] <- Varw + tcrossprod(w[, j])
        logdetA[j] <- determinant(A)$mod
      } else {
        # Nystrom approximation
        K.mm <- Hlam.matsq[[j]][1:Nystrom$m, 1:Nystrom$m]
        eigenK.mm <- iprior::eigenCpp(K.mm)
        V <- Hlam.matsq[[j]][, 1:Nystrom$m] %*% eigenK.mm$vec
        u <- eigenK.mm$val
        u.Vt <- t(V) * u
        D <- u.Vt %*% V + diag(1, Nystrom$m)
        E <- solve(D, u.Vt, tol = 1e-18)
        # see https://stackoverflow.com/questions/22134398/mahalonobis-distance-in-r-error-system-is-computationally-singular
        # see https://stackoverflow.com/questions/21451664/system-is-computationally-singular-error
        w[, j] <- as.numeric(a - V %*% (E %*% a))
        W[[j]] <- (diag(1, n) - V %*% E) + tcrossprod(w[, j])
        logdetA[j] <- determinant(A)$mod
      }
    }
    if (isTRUE(common.intercept) && isTRUE(common.RKHS.scale)) {
      w <- matrix(w[, 1], nrow = n, ncol = m)
      W <- rep(list(W[[1]]), m)
      logdetA <- rep(logdetA[1], m)
    }

    # Update lambda ------------------------------------------------------------
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

    # Update H.lam and H.lam.sq ------------------------------------------------
    lambdaExpand_mult(env = iprobit.env)
    HlamFn_mult(env = iprobit.env)
    HlamsqFn_mult(env = iprobit.env)

    # Update alpha -------------------------------------------------------------
    alpha <- apply(ystar - mapply("%*%", Hlam.mat, split(w, col(w))), 2, mean)
    if (isTRUE(common.intercept)) alpha <- rep(mean(alpha), m)

    # Calculate lower bound ----------------------------------------------------
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
    lb[niter + 1] <- lb.ystar + lb.w + lb.lambda + lb.alpha

    # Calculate fitted values and error rate -----------------------------------
    ystar <- rep(alpha, each = n) + mapply("%*%", Hlam.mat, split(w, col(w)))
    fitted.values <- predict_iprobit_mult(y, y.levels, ystar)
    error.rates[niter + 1] <- fitted.values$train.error
    brier.scores[niter + 1] <- fitted.values$brier.score

    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, niter)
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

  # Clean up and close ---------------------------------------------------------
  lambda <- matrix(lambda[1:l, ], ncol = m, nrow = l)
  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar, w = w, lambda = lambda, alpha = alpha,
              lower.bound = lb[!is.na(lb)], ipriorKernel = NULL,
              se.alpha = se.alpha, se.lambda = se.lambda, se.ystar = se.ystar,
              y.levels = y.levels, start.time = start.time, end.time = end.time,
              time = time.taken, stop.crit = stop.crit, niter = niter,
              maxit = maxit, fitted.values = fitted.values,
              error = error.rates[!is.na(error.rates)],
              brier = brier.scores[!is.na(brier.scores)])
  class(res) <- c("iprobitMod", "iprobitMod_mult")
  res
}
