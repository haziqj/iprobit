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

iprobit_mult <- function(mod, maxit = 10, stop.crit = 1e-5, silent = FALSE,
                         alpha0 = NULL, theta0 = NULL, w0 = NULL,
                         common.intercept = FALSE, common.RKHS.scale = FALSE) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(mod, iprobit.env)
  list2env(BlockBStuff, iprobit.env)
  environment(BlockB) <- iprobit.env
  environment(get_Hlamsq) <- iprobit.env
  environment(loop_logical) <- iprobit.env
  y <- as.numeric(factor(y))  # as.factor then as.numeric to get y = 1, 2, ...
  m <- length(y.levels)
  nm <- n * m
  maxit <- max(1, maxit)

  # Initialise -----------------------------------------------------------------
  if (is.null(alpha0)) {
    if (isTRUE(common.intercept)) alpha0 <- rep(rnorm(1), m)
    else alpha0 <- rnorm(m)
  }
  if (is.null(theta0)) {
    if (isTRUE(common.RKHS.scale)) theta0 <- rep(rnorm(p), m)
    else theta0 <- rnorm(p * m)
  }
  if (p == 1) lambda0 <- exp(theta0)
  else lambda0 <- theta0
  if (is.null(w0)) w0 <- matrix(0, ncol = m, nrow = n)

  alpha <- rep(1, m); alpha[] <- alpha0  # sometimes it is convenient to set alpha0 = 1
  theta <- lambda <- ct <- dt <- matrix(lambda0, ncol = m, nrow = p)
  lambdasq <- lambda ^ 2
  Hl <- iprior::.expand_Hl_and_lambda(Hl, rep(1, p), intr, intr.3plus)$Hl  # expand Hl
  lambda <- expand_lambda(lambda[1:p, , drop = FALSE], intr, intr.3plus)
  lambdasq <- expand_lambda(lambdasq[1:p, , drop = FALSE], intr, intr.3plus)
  Hlam   <- get_Hlam(mod, lambda[1:p, , drop = FALSE], theta.is.lambda = TRUE)
  Hlamsq <- get_Hlamsq()
  w <- f.tmp <- ystar <- w0
  Varw <- W <- list(NULL)
  f.tmp <- rep(alpha, each = n) + mapply("%*%", Hlam, split(w, col(w)))
  logdetA <- rep(NA, m)
  w.ind <- seq_len(
    ifelse(isTRUE(common.intercept) && isTRUE(common.RKHS.scale), 1, m)
  )
  niter <- 0
  lb <- train.error <- train.brier <- test.error <- test.brier <- rep(NA, maxit)
  logClb <- rep(NA, n)

  # The variational EM loop ----------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (loop_logical()) {  # see loop_logical() function in iprobit_helper.R
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
      A <- Hlamsq[[j]] + diag(1, n)
      a <- as.numeric(crossprod(Hlam[[j]], ystar[, j] - alpha[j]))
      eigenA <- iprior::eigenCpp(A)
      V <- eigenA$vec
      u <- eigenA$val + 1e-8  # ensure positive eigenvalues
      uinv.Vt <- t(V) / u
      w[, j] <- as.numeric(crossprod(a, V) %*% uinv.Vt)
      Varw[[j]] <- iprior::fastVDiag(V, 1 / u)  # V %*% uinv.Vt
      W[[j]] <- Varw[[j]] + tcrossprod(w[, j])
      logdetA[j] <- determinant(A)$mod
    }
    if (isTRUE(common.intercept) && isTRUE(common.RKHS.scale)) {
      w <- matrix(w[, 1], nrow = n, ncol = m)
      W <- rep(list(W[[1]]), m)
      logdetA <- rep(logdetA[1], m)
    }

    # Update lambda ------------------------------------------------------------
    for (k in 1:p) {
      for (j in 1:m) {
        lambda   <- expand_lambda(lambda[1:p, , drop = FALSE], intr)
        lambdasq <- expand_lambda(lambdasq[1:p, , drop = FALSE], intr)
        BlockB(k, lambda[, j])
        ct[k, j] <- sum(Psql[[k]] * W[[j]])
        dt[k, j] <- as.numeric(
          crossprod(ystar[, j] - alpha[j], Pl[[k]]) %*% w[, j] -
            sum(Sl[[k]] * W[[j]]) / 2
        )
      }
      if (isTRUE(common.RKHS.scale)) {
        lambda[k, ] <- rep(sum(dt[k, ]) / sum(ct[k, ]), m)
        lambdasq[k, ] <- rep(1 / sum(ct[k, ]) + lambda[k, 1] ^ 2, m)
      } else {
        lambda[k, ] <- dt[k, ] / ct[k, ]
        lambdasq[k, ] <- 1 / ct[k, ] + lambda[k, ] ^ 2
      }
    }

    # Update H.lam and H.lam.sq ------------------------------------------------
    lambda   <- expand_lambda(lambda[1:p, , drop = FALSE], intr)
    lambdasq <- expand_lambda(lambdasq[1:p, , drop = FALSE], intr)
    Hlam     <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)
    Hlamsq   <- get_Hlamsq()

    # Update alpha -------------------------------------------------------------
    alpha <- apply(ystar - mapply("%*%", Hlam, split(w, col(w))), 2, mean)
    if (isTRUE(common.intercept)) alpha <- rep(mean(alpha), m)

    # theta --------------------------------------------------------------------
    if (p == 1) theta <- log(lambda[1, , drop = FALSE])
    else theta <- lambda[1:p, , drop = FALSE]

    # Calculate lower bound ----------------------------------------------------
    lb.ystar <- sum(logClb)
    lb.w <- 0.5 * (nm - sum(sapply(W, function(x) sum(diag(x)))) - sum(logdetA))
    if (isTRUE(common.RKHS.scale))
      lb.lambda <- (p / 2) * (1 + log(2 * pi)) - sum(log(apply(ct, 1, sum))) / 2
    else
      lb.lambda <- (m / 2) * (p * (1 + log(2 * pi)) - sum(log(ct)) / m)
    if (isTRUE(common.intercept))
      lb.alpha <- 0.5 * (1 + log(2 * pi) - log(nm))
    else
      lb.alpha <- (m / 2) * (1 + log(2 * pi) - log(n))
    lb[niter + 1] <- lb.ystar + lb.w + lb.lambda + lb.alpha

    # Calculate fitted values and error rate -----------------------------------
    f.tmp <- rep(alpha, each = n) + mapply("%*%", Hlam, split(w, col(w)))
    fitted.values <- probs_yhat_error(y, y.levels, f.tmp)
    train.error[niter + 1] <- fitted.values$error
    train.brier[niter + 1] <- fitted.values$brier
    fitted.test <- NULL
    if (iprior::.is.ipriorKernel_cv(mod)) {
      ystar.test <- calc_ystar(mod, mod$Xl.test, alpha, theta, w)
      fitted.test <- probs_yhat_error(y.test, y.levels, ystar.test)
      test.error[niter + 1] <- fitted.test$error
      test.brier[niter + 1] <- fitted.test$brier
    }

    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, niter)
  }

  end.time <- Sys.time()
  time.taken <- iprior::as.time(end.time - start.time)

  # Calculate posterior s.d. and prepare param table ---------------------------
  param.full <- theta_to_param.full(theta, alpha, mod)
  if (isTRUE(common.intercept)) {
    alpha <- unique(alpha)  # since they're all the same
    se.alpha <- sqrt(1 / nm)
  } else {
    se.alpha <- rep(sqrt(1 / n), m)
  }
  lambda <- lambda[1:p, , drop = FALSE]
  if (isTRUE(common.RKHS.scale)) {
    lambda <- apply(lambda, 1, unique)  # since they're all the same
    se.lambda <- matrix(sqrt(1 / apply(ct, 1, sum)), ncol = m, nrow = p)[, 1]
  } else {
    lambda <- c(t(lambda))
    se.lambda <- c(t(matrix(sqrt(1 / ct[1:p, ]), ncol = m, nrow = p)))
  }
  param.summ <- data.frame(
    Mean    = c(alpha, lambda),
    S.D.    = c(se.alpha, se.lambda),
    `2.5%`  = c(alpha, lambda) - qnorm(0.975) * c(se.alpha, se.lambda),
    `97.5%` = c(alpha, lambda) + qnorm(0.975) * c(se.alpha, se.lambda)
  )
  colnames(param.summ) <- c("Mean", "S.D.", "2.5%", "97.5%")
  rownames(param.summ) <- c(get_names(mod, "intercept", !common.intercept),
                            get_names(mod, "lambda", !common.RKHS.scale))

  # Clean up and close ---------------------------------------------------------
  convergence <- niter == maxit
  if (!silent) {
    close(pb)
    if (convergence) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  list(theta = theta, param.full = param.full, param.summ = param.summ, w = w,
       Varw = Varw, lower.bound = as.numeric(na.omit(lb)), niter = niter,
       start.time = start.time, end.time = end.time, time = time.taken,
       fitted.values = fitted.values, test = fitted.test,
       train.error = as.numeric(na.omit(train.error)),
       train.brier = as.numeric(na.omit(train.brier)),
       test.error = as.numeric(na.omit(test.error)),
       test.brier = as.numeric(na.omit(test.brier)), convergence = convergence)
}

# else {
#   # Nystrom approximation
#   K.mm <- Hlamsq[[j]][1:Nystrom$m, 1:Nystrom$m]
#   eigenK.mm <- iprior::eigenCpp(K.mm)
#   V <- Hlamsq[[j]][, 1:Nystrom$m] %*% eigenK.mm$vec
#   u <- eigenK.mm$val
#   u.Vt <- t(V) * u
#   D <- u.Vt %*% V + diag(1, Nystrom$m)
#   E <- solve(D, u.Vt, tol = 1e-18)
#   # see https://stackoverflow.com/questions/22134398/mahalonobis-distance-in-r-error-system-is-computationally-singular
#   # see https://stackoverflow.com/questions/21451664/system-is-computationally-singular-error
#   w[, j] <- as.numeric(a - V %*% (E %*% a))
#   W[[j]] <- (diag(1, n) - V %*% E) + tcrossprod(w[, j])
#   logdetA[j] <- determinant(A)$mod
# }
