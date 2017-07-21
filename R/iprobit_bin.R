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
iprobit_bin <- function(ipriorKernel, maxit = 100, stop.crit = 1e-5,
                        silent = FALSE, alpha0 = NULL, lambda0 = NULL,
                        w0 = NULL) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(ipriorKernel, iprobit.env)
  list2env(BlockBstuff, iprobit.env)
  list2env(model, iprobit.env)
  environment(BlockB) <- iprobit.env
  environment(lambdaExpand_bin) <- iprobit.env
  environment(HlamFn) <- iprobit.env
  environment(HlamsqFn) <- iprobit.env
  y <- Y
  maxit <- max(1, maxit)

  # Initialise -----------------------------------------------------------------
  if (is.null(alpha0)) alpha0 <- rnorm(1)
  if (is.null(lambda0)) lambda0 <- abs(rnorm(l))
  if (is.null(w0)) w0 <- rep(0, n)
  lambda <- ct <- lambda0
  lambda.sq <- lambda0 ^ 2
  lambdaExpand_bin(env = iprobit.env)
  HlamFn(env = iprobit.env)
  HlamsqFn(env = iprobit.env)
  alpha <- alpha0
  w <- w0

  # Variational lower bound and loopy stuff ------------------------------------
  niter <- 0
  lb <- rep(NA, maxit)
  lb.const <- (n + 1 + l - log(n) + (l + 1) * log(2 * pi)) / 2
  loop.logical <- function() {
    lb.diff <- (lb[niter] - lb[niter - 1])
    ifelse(length(lb.diff) == 0, TRUE,
           ifelse(is.na(lb.diff), niter != maxit,
                  (niter != maxit) && (lb.diff > stop.crit)))
  }

  # The variational EM loop ----------------------------------------------------
  if (maxit == 1) silent <- TRUE
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 1)
  start.time <- Sys.time()

  while (loop.logical()) {
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
    if (!isNystrom(ipriorKernel)) {
      eigenA <- iprior::eigenCpp(A)
      V <- eigenA$vec
      u <- eigenA$val + 1e-8  # ensure positive eigenvalues
      uinv.Vt <- t(V) / u
      w <- as.numeric(crossprod(a, V) %*% uinv.Vt)
      Varw <- iprior::fastVDiag(V, 1 / u)  # V %*% uinv.Vt
      W <- Varw + tcrossprod(w)
    } else {
      # Nystrom approximation
      K.mm <- Hlam.matsq[1:Nystrom$m, 1:Nystrom$m]
      eigenK.mm <- eigen(K.mm)
      V <- Hlam.matsq[, 1:Nystrom$m] %*% eigenK.mm$vec
      u <- eigenK.mm$val
      u.Vt <- t(V) * u
      D <- u.Vt %*% V + diag(1, Nystrom$m)
      E <- solve(D, u.Vt)
      w <- as.numeric(a - V %*% (E %*% a))
      W <- (diag(1, n) - V %*% E) + tcrossprod(w)
    }

    # Update lambda ------------------------------------------------------------
    for (k in 1:l) {
      lambdaExpand_bin(env = iprobit.env)
      BlockB(k)
      ct[k] <- sum(Psql[[k]] * W)
      dt <- as.numeric(
        crossprod(ystar - alpha, Pl[[k]]) %*% w - sum(Sl[[k]] * W) / 2
      )
      lambda[k] <- dt / ct[k]
      lambda.sq[k] <- 1 / ct[k] + (dt / ct[k]) ^ 2
    }
    lambdaExpand_bin(env = iprobit.env)
    HlamFn(env = iprobit.env)
    HlamsqFn(env = iprobit.env)

    # Update alpha -------------------------------------------------------------
    alpha <- mean(ystar - Hlam.mat %*% w)

    # Calculate lower bound ----------------------------------------------------
    lb[niter + 1] <- lb.const +
      sum(pnorm(eta[y == 2], log.p = TRUE)) +
      sum(pnorm(-eta[y == 1], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus + sum(log(ct))) / 2

    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, niter)
  }

  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Calculate standard errors from posterior variance --------------------------
  se.alpha <- sqrt(1 / n)
  se.lambda <- sqrt(1 / ct[1:l])
  se.ystar <- NA #iprobitSE(y = y, eta = eta, thing1 = thing1, thing0 = thing0)

  # Clean up and close ---------------------------------------------------------
  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar, w = w, lambda = lambda[1:l], alpha = alpha,
              lower.bound = lb, ipriorKernel = ipriorKernel,
              se.alpha = se.alpha, se.lambda = se.lambda, se.ystar = se.ystar,
              y.levels = y.levels, start.time = start.time,
              end.time = end.time, time = time.taken,
              stop.crit = stop.crit, niter = niter, maxit = maxit)
  class(res) <- c("iprobitMod", "iprobitMod_bin")
  res
}
