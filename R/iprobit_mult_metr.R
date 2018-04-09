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

# theta.samp would now be a list of length n.samp and each component of this
# list is the matrix theta.samp.
#
# Similarly Hlaml and Hlamsql are lists of length n.samp, and each component of
# the list is another list of length  m (no. of classes) where the values
# Hlam_theta for each theta.samp are stored.

iprobit_mult_metr <- function(mod, maxit = 5, stop.crit = 1e-5, silent = FALSE,
                              alpha0 = NULL, theta0 = NULL, w0 = NULL,
                              n.samp = 100, sd.samp = 0.15, thin.samp = 1,
                              seed = NULL) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(mod, iprobit.env)
  environment(loop_logical) <- iprobit.env
  y <- as.numeric(factor(y))  # as.factor then as.numeric to get y = 1, 2, ...
  m <- length(y.levels)
  nm <- n * m
  maxit <- max(1, maxit)
  thin.seq <- seq_len(n.samp)[seq_len(n.samp) %% thin.samp == 0]

  # Initialise -----------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  if (is.null(alpha0)) alpha0 <- rnorm(m)
  if (is.null(theta0)) theta0 <- rep(abs(rnorm(p)), m)
  if (is.null(w0)) w0 <- matrix(0, ncol = m, nrow = n)

  alpha <- alpha0
  theta.samp <- Hlaml <- Hlamsql <- list(NULL)
  theta <- matrix(theta0, ncol = m, nrow = thetal$n.theta)
  Hlaml[[1]] <- Hlam <- get_Hlam(mod, theta)
  Hlamsql[[1]] <- Hlamsq <- lapply(Hlaml[[1]], iprior::fastSquare)
  theta.samp[[1]] <- theta

  w <- f.tmp <- ystar <- w0
  Varw <- W <- list(NULL)
  f.tmp <- rep(alpha, each = n) + mapply("%*%", Hlam, split(w, col(w)))
  logdetA <- rep(NA, m)
  niter <- 0
  lb <- train.error <- train.brier <- test.error <- test.brier <- rep(NA, maxit)
  logClb <- rep(NA, n)
  log.q.star <- rep(NA, m)
  log.q <- matrix(NA, ncol = m, nrow = n.samp)
  log.q[1, ] <- sapply(Hlamsq, function(x) -0.5 * sum(diag(x)))
  acc.rate <- matrix(NA, nrow = maxit, ncol = m)

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
    for (j in 1:m) {
      A <- Hlamsq + diag(1, n)
      a <- as.numeric(crossprod(Hlam, ystar[, j] - alpha[j]))
      eigenA <- iprior::eigenCpp(A)
      V <- eigenA$vec
      u <- eigenA$val + 1e-8  # ensure positive eigenvalues
      uinv.Vt <- t(V) / u
      w[, j] <- as.numeric(crossprod(a, V) %*% uinv.Vt)
      Varw[[j]] <- iprior::fastVDiag(V, 1 / u)  # V %*% uinv.Vt
      W[[j]] <- Varw[[j]] + tcrossprod(w[, j])
      logdetA[j] <- determinant(A)$mod
    }

    # Update theta -------------------------------------------------------------
    count <- rep(0, m)
    for (t in seq_len(n.samp - 1)) {
      theta.star <- theta.current <- theta.samp[[t]]
      theta.star[] <- rnorm(thetal$n.theta, mean = theta.current[, 1], sd = sd.samp)
      Hlam.star <- get_Hlam(mod, theta.star)
      Hlamsq.star <- iprior::fastSquare(Hlam.star)
      log.q.star[] <- as.numeric(-0.5 * (
        sum(Hlamsq.star * W[[1]]) -
          2 * crossprod(ystar[, 1] - alpha[1], Hlam.star) %*% w[, 1]
      ))
      log.prob.acc <- log.q.star - log.q[t, ]
      prob.acc <- exp(log.prob.acc)
      prob.acc[prob.acc > 1] <- 1
      if (runif(1) < prob.acc[j]) {
        theta.samp[[t + 1]] <- theta.star
        Hlaml[[t + 1]] <- Hlam.star
        Hlamsql[[t + 1]] <- Hlamsq.star
        log.q[t + 1, ] <- log.q.star
        count <- count + 1
      } else {
        theta.samp[[t + 1]] <- theta.current
        Hlaml[[t + 1]] <- Hlaml[[t]]
        Hlamsql[[t + 1]] <- Hlamsql[[t]]
        log.q[t + 1, ] <- log.q[t, ]
      }
    }
    acc.rate[niter, ] <- count / n.samp
    theta <- Reduce("+", theta.samp[thin.seq]) / length(theta.samp[thin.seq])

    # Update Hlam and Hlamsq ---------------------------------------------------
    for (j in 1:m) {
      tmp <- lapply(Hlaml, function(x) x[[j]])
      Hlam[[j]] <- Reduce("+", tmp[thin.seq]) / length(tmp[thin.seq])
      tmp <- lapply(Hlamsql, function(x) x[[j]])
      Hlamsq[[j]] <- Reduce("+", tmp[thin.seq]) / length(tmp[thin.seq])
    }

    # Update alpha -------------------------------------------------------------
    alpha <- apply(ystar - mapply("%*%", Hlam, split(w, col(w))), 2, mean)
    alpha <- alpha - mean(alpha)
    # if (isTRUE(common.intercept)) alpha <- rep(mean(alpha), m)

    # Calculate lower bound ----------------------------------------------------
    lb.ystar <- sum(logClb)
    lb.w <- 0.5 * (nm - sum(sapply(W, function(x) sum(diag(x)))) - sum(logdetA))
    lb.lambda <- -sum(apply(log.q, 2, mean))
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
  se.alpha <- rep(sqrt(1 / n), m)
  alpha.summ <- cbind(alpha, se.alpha, alpha - qnorm(0.975) * se.alpha,
                      alpha + qnorm(0.975) * se.alpha)
  theta.summ <- theta.samp_to_summ_mult(theta.samp, mod)
  colnames(alpha.summ) <- colnames(theta.summ)
  param.summ <- rbind(alpha.summ, theta.summ)
  rownames(param.summ) <- c(get_names(mod, "intercept", TRUE),
                            get_names(mod, c("lambda", "hurst", "lengthscale",
                                             "offset"), FALSE))

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
       test.brier = as.numeric(na.omit(test.brier)), convergence = convergence,
       acc.rate = as.numeric(na.omit(acc.rate)),
       theta.samp = theta.samp[thin.seq])
}

theta.samp_to_summ_mult <- function(theta.samp, object) {
  m <- ncol(theta.samp[[1]])
  theta.samp_col <- function(x, j) x[, j]
  res <- NULL
  for (j in seq_len(m)) {
    tmp <- t(sapply(theta.samp, theta.samp_col, j = j))
    res[[j]] <- theta.samp_to_summ(tmp, object)
  }
  res.lambda <- do.call(rbind, lapply(res, function(x) x[grep("lambda", rownames(x)), ]))
  res.hurst <- do.call(rbind, lapply(res, function(x) x[grep("hurst", rownames(x)), ]))
  res.lengthscale <- do.call(rbind, lapply(res, function(x) x[grep("lengthscale", rownames(x)), ]))
  res.offset <- do.call(rbind, lapply(res, function(x) x[grep("offset", rownames(x)), ]))

  res <- rbind(
    res.lambda, res.hurst, res.lengthscale, res.offset
  )
  unq.ind <- match(unique(res[, 1]), res[, 1])
  # rownames(res) <- get_names(mod, c("lambda", "hurst", "lengthscale", "offset"))
  res[unq.ind, ]
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
