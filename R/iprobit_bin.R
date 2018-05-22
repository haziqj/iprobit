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
#   GNU General Public License for more df.tmpils.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

iprobit_bin <- function(mod, maxit = 100, stop.crit = 1e-5, silent = FALSE,
                        alpha0 = NULL, theta0 = NULL, w0 = NULL, w.only = FALSE,
                        int.only = FALSE) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(mod, iprobit.env)
  list2env(BlockBStuff, iprobit.env)
  environment(BlockB) <- iprobit.env
  environment(get_Hlamsq) <- iprobit.env
  environment(loop_logical) <- iprobit.env
  y <- as.numeric(factor(y))  # as.factor then as.numeric to get y = 1, 2, ...
  maxit <- max(1, maxit)  # cannot have maxit <= 0

  # Initialise -----------------------------------------------------------------
  if (is.null(alpha0)) alpha0 <- rnorm(1)
  if (is.null(theta0)) theta0 <- rnorm(p)
  if (is.null(w0)) w0 <- rep(0, n)
  if (p == 1) lambda0 <- exp(theta0)  # ensures positive values
  else lambda0 <- theta0

  alpha <- alpha0
  lambda <- ct <- lambda0
  lambdasq <- lambda0 ^ 2
  Hl <- iprior::.expand_Hl_and_lambda(Hl, rep(1, p), intr, intr.3plus)$Hl  # expand Hl
  lambda <- expand_lambda(lambda[1:p], intr, intr.3plus)
  lambdasq <- expand_lambda(lambdasq[1:p], intr, intr.3plus)
  Hlam   <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)
  w <- w0
  Varw <- NA
  f.tmp <- as.numeric(alpha + Hlam %*% w)
  niter <- 0
  lb <- train.error <- train.brier <- test.error <- test.brier <- rep(NA, maxit)

  # The variational EM loop ----------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (loop_logical()) {  # see loop_logical() function in iprobit_helper.R
    # Update ystar -------------------------------------------------------------
    thing <- rep(NA, n)
    thing1 <- exp(  # phi(f.tmp) / Phi(f.tmp)
      dnorm(f.tmp[y == 2], log = TRUE) - pnorm(f.tmp[y == 2], log.p = TRUE)
    )
    thing0 <- -exp(  # -1 * {phi(f.tmp) / Phi(-f.tmp)}
      dnorm(f.tmp[y == 1], log = TRUE) - pnorm(-f.tmp[y == 1], log.p = TRUE)
    )
    thing[y == 2] <- thing1
    thing[y == 1] <- thing0
    ystar <- f.tmp + thing

    # Update w -----------------------------------------------------------------
    iprior::.eigen_Hlam(Hlam, environment())  # assign u and V to environment
    z <- u ^ 2 + 1  # eigenvalues of Vy
    zinv.Vt <- t(V) / z
    Vy.inv.y <- as.numeric(crossprod(ystar - alpha, V) %*% zinv.Vt)
    w <- Hlam %*% Vy.inv.y
    W <- V %*% zinv.Vt + tcrossprod(w)
    Varw <- V %*% zinv.Vt

    if (!isTRUE(w.only)) {
      # Update lambda ----------------------------------------------------------
      for (k in 1:p) {
        lambda   <- expand_lambda(lambda[1:p], intr)
        BlockB(k)  # Updates Pl, Psql, and Sl in environment
        ct[k] <- sum(Psql[[k]] * W)
        dt <- as.numeric(
          crossprod(ystar - alpha, Pl[[k]]) %*% w - sum(Sl[[k]] * W) / 2
        )
        if (!isTRUE(int.only)) {
          lambda[k] <- dt / ct[k]
        }
      }
      lambda   <- expand_lambda(lambda[1:p], intr, intr.3plus)
      Hlam     <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)

      # Update alpha -----------------------------------------------------------
      alpha <- mean(ystar - Hlam %*% w)
    }

    # Calculate lower bound ----------------------------------------------------
    lb[niter + 1] <- n / 2  +
      sum(pnorm(f.tmp[y == 2], log.p = TRUE)) +
      sum(pnorm(-f.tmp[y == 1], log.p = TRUE)) -
      0.5 * sum(diag(W)) +
      0.5 * as.numeric(determinant(Varw)$mod)

    # theta --------------------------------------------------------------------
    if (p == 1) theta <- log(lambda[1])
    else theta <- lambda[1:p]
    theta <- matrix(theta, ncol = 2, nrow = length(theta))

    # Calculate fitted values and error rates ----------------------------------
    f.tmp <- as.numeric(alpha + Hlam %*% w)  # E[ystar|y,theta]
    f.var.tmp <- diag(Hlam %*% Varw %*% Hlam) + 1  # diag(Var[ystar|y,theta])
    fitted.values <- probs_yhat_error(y, y.levels, f.tmp / sqrt(f.var.tmp))
    train.error[niter + 1] <- fitted.values$error
    train.brier[niter + 1] <- fitted.values$brier
    fitted.test <- NULL
    if (iprior::.is.ipriorKernel_cv(mod)) {
      ystar.test <- calc_ystar(mod, mod$Xl.test, alpha, theta, w, Varw = Varw)
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
  param.summ <- data.frame(
    Mean    = c(alpha, lambda[1:p]),
    S.D.    = sqrt(rep(0, p + 1)),
    `2.5%`  = sqrt(rep(0, p + 1)),
    `97.5%` = sqrt(rep(0, p + 1))
  )
  colnames(param.summ) <- c("Mean", "S.D.", "2.5%", "97.5%")
  rownames(param.summ) <- get_names(mod, c("intercept", "lambda"), FALSE)

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
