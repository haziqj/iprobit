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

#' @export
iprobit_bin <- function(mod, maxit = 100, stop.crit = 1e-5, silent = FALSE,
                        alpha0 = NULL, lambda0 = NULL, w0 = NULL) {
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
  if (is.null(lambda0)) lambda0 <- abs(rnorm(p))
  if (is.null(w0)) w0 <- rep(0, n)
  lambda <- ct <- lambda0
  lambdasq <- lambda0 ^ 2
  Hl <- iprior::.expand_Hl_and_lambda(Hl, rep(1, p), intr, intr.3plus)$Hl  # expand Hl
  lambda <- expand_lambda(lambda[1:p], intr, intr.3plus)
  lambdasq <- expand_lambda(lambdasq[1:p], intr, intr.3plus)
  Hlam   <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)
  Hlamsq <- get_Hlamsq()
  alpha <- alpha0
  w <- w0
  Varw <- NA
  f.tmp <- as.numeric(alpha + Hlam %*% w)
  niter <- 0
  lb <- train.error <- train.brier <- test.error <- test.brier <- rep(NA, maxit)
  lb.const <- (n + 1 + p - log(n) + (p + 1) * log(2 * pi)) / 2

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
    A <- Hlamsq + diag(1, n)
    a <- as.numeric(crossprod(Hlam, ystar - alpha))
    eigenA <- iprior::eigenCpp(A)
    V <- eigenA$vec
    u <- eigenA$val
    uinv.Vt <- t(V) / u
    w <- as.numeric(crossprod(a, V) %*% uinv.Vt)
    Varw <- iprior::fastVDiag(V, 1 / u)  # V %*% uinv.Vt
    W <- Varw + tcrossprod(w)

    # Update lambda ------------------------------------------------------------
    for (k in 1:p) {
      lambda   <- expand_lambda(lambda[1:p], intr)
      lambdasq <- expand_lambda(lambdasq[1:p], intr)
      BlockB(k)  # Updates Pl, Psql, and Sl in environment
      ct[k] <- sum(Psql[[k]] * W)
      dt <- as.numeric(
        crossprod(ystar - alpha, Pl[[k]]) %*% w - sum(Sl[[k]] * W) / 2
      )
      lambda[k] <- dt / ct[k]
      lambdasq[k] <- 1 / ct[k] + (dt / ct[k]) ^ 2
    }
    lambda   <- expand_lambda(lambda[1:p], intr, intr.3plus)
    lambdasq <- expand_lambda(lambdasq[1:p], intr, intr.3plus)
    Hlam     <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)
    Hlamsq   <- get_Hlamsq()

    # Update alpha -------------------------------------------------------------
    alpha <- mean(ystar - Hlam %*% w)

    # Calculate lower bound ----------------------------------------------------
    lb[niter + 1] <- lb.const +
      sum(pnorm(f.tmp[y == 2], log.p = TRUE)) +
      sum(pnorm(-f.tmp[y == 1], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus + sum(log(ct))) / 2

    # theta --------------------------------------------------------------------
    theta <- hyperparam_to_theta(lambda = lambda[1:p])

    # Calculate fitted values and error rates ----------------------------------
    f.tmp <- as.numeric(alpha + Hlam %*% w)
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
  hyperparam <- c(alpha, lambda[1:p])
  se <- sqrt(c(1 / n, 1 / ct[1:p]))  # c(alpha, lambda)
  param.summ <- data.frame(
    Mean    = hyperparam,
    S.D.    = se,
    `2.5%`  = hyperparam - qnorm(0.975) * se,
    `97.5%` = hyperparam + qnorm(0.975) * se
  )
  colnames(param.summ) <- c("Mean", "S.D.", "2.5%", "97.5%")
  if (length(lambda) > 1) {
    rownames(param.summ) <- c("alpha", paste0("lambda[", seq_along(lambda), "]"))
  } else {
    rownames(param.summ) <- c("alpha", "lambda")
  }

  # Clean up and close ---------------------------------------------------------
  convergence <- niter == maxit
  if (!silent) {
    close(pb)
    if (convergence) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  list(theta = theta, param.summ = param.summ, w = w, Varw = Varw,
       lower.bound = as.numeric(na.omit(lb)), niter = niter,
       start.time = start.time, end.time = end.time, time = time.taken,
       fitted.values = fitted.values, test = fitted.test,
       train.error = as.numeric(na.omit(train.error)),
       train.brier = as.numeric(na.omit(train.brier)),
       test.error = as.numeric(na.omit(test.error)),
       test.brier = as.numeric(na.omit(test.brier)), convergence = convergence)
}

# if (!isNystrom(mod)) {
# } else {
#   # Nystrom approximation
#   K.mm <- Hlamsq[1:Nystrom$m, 1:Nystrom$m]
#   eigenK.mm <- iprior::eigenCpp(K.mm)
#   V <- Hlamsq[, 1:Nystrom$m] %*% eigenK.mm$vec
#   u <- eigenK.mm$val
#   u.Vt <- t(V) * u
#   D <- u.Vt %*% V + diag(1, Nystrom$m)
#   E <- solve(D, u.Vt, tol = 1e-18)
#   # see https://stackoverflow.com/questions/22134398/mahalonobis-distance-in-r-error-system-is-computationally-singular
#   # see https://stackoverflow.com/questions/21451664/system-is-computationally-singular-error
#   w <- as.numeric(a - V %*% (E %*% a))
#   W <- (diag(1, n) - V %*% E) + tcrossprod(w)
# }
