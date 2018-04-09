iprobit_bin_is <- function(mod, maxit = 20, stop.crit = 1e-5, silent = FALSE,
                             alpha0 = NULL, theta0 = NULL, w0 = NULL,
                             n.samp = 100, sd.samp = 1) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  mod$estl$est.psi <- FALSE
  list2env(mod, iprobit.env)
  environment(loop_logical) <- iprobit.env
  maxit <- max(1, maxit)  # cannot have maxit <= 0
  y <- as.numeric(factor(y))  # as.factor then as.numeric to get y = 1, 2, ...

  # Initialise -----------------------------------------------------------------
  if (is.null(alpha0)) alpha0 <- rnorm(1)
  if (is.null(theta0)) theta0 <- rnorm(thetal$n.theta)
  if (is.null(w0)) w0 <- rep(0, n)
  alpha <- alpha0
  theta <- theta0
  Hlam <- iprior::.get_Hlam(mod, theta)
  Hlamsq <- iprior::fastSquare(Hlam)
  theta.samp <- matrix(NA, nrow = n.samp, ncol = thetal$n.theta)
  theta.samp[1, ] <- theta
  w <- w0
  Varw <- NA
  eta <- as.numeric(alpha + Hlam %*% w)
  niter <- 0
  acc.rate <- lb <- train.error <- train.brier <- test.error <- test.brier <- rep(NA, maxit)
  lb.const <- (n + 1 + p - log(n) + (p + 1) * log(2 * pi)) / 2

  # The variational EM loop ----------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (niter != maxit) {  # see loop_logical() function in iprobit_helper.R
    # Update ystar -------------------------------------------------------------
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
    A <- Hlamsq + diag(1, n)
    a <- as.numeric(crossprod(Hlam, ystar - alpha))
    eigenA <- iprior::eigenCpp(A)
    V <- eigenA$vec
    u <- eigenA$val
    uinv.Vt <- t(V) / u
    w <- as.numeric(crossprod(a, V) %*% uinv.Vt)
    Varw <- iprior::fastVDiag(V, 1 / u)  # V %*% uinv.Vt
    W <- Varw + tcrossprod(w)

    # Update theta -------------------------------------------------------------
    theta.samp <- mvtnorm::rmvnorm(n.samp, mean = theta,
                                   sigma = diag(sd.samp ^ 2, length(theta)))
    tmp <- split(theta.samp, rep(seq_len(nrow(theta.samp)),
                                 each = ncol(theta.samp)))
    weights.isl <- lapply(tmp, function(x) {
      res <- log_q_eta(x, mod, W, ystar - alpha, w) -
        mvtnorm::dmvnorm(x, mean = theta, sigma = diag(sd.samp ^ 2, length(x)),
                         log = TRUE)
      exp(res)
    })
    Hlaml <- lapply(tmp, iprior::.get_Hlam, object = mod)
    Hlamsql <- lapply(tmp, function(x) {
      iprior::fastSquare(iprior::.get_Hlam(x, object = mod))
    })
    weights.is <- Reduce("+", weights.isl)
    theta <- apply(theta.samp, 2, mean) / weights.is
    Hlam <- Reduce("+", Hlaml) / weights.is
    Hlamsq <- Reduce("+", Hlamsql) / weights.is
    # Hlam <- iprior::.get_Hlam(mod, theta)
    # Hlamsq <- iprior::fastSquare(Hlam)

    # Update alpha -------------------------------------------------------------
    alpha <- mean(ystar - Hlam %*% w)

    # Calculate lower bound ----------------------------------------------------
    lb[niter + 1] <- lb.const +
      sum(pnorm(eta[y == 2], log.p = TRUE)) +
      sum(pnorm(-eta[y == 1], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus ) / 2

    # Calculate fitted values and error rates ----------------------------------
    eta <- as.numeric(alpha + Hlam %*% w)
    fitted.values <- probs_yhat_error(y, y.levels, eta)
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
  # hyperparam <- c(alpha, lambda[1:p])
  # se <- sqrt(c(1 / n, 1 / ct[1:p]))  # c(alpha, lambda)
  # param.summ <- data.frame(
  #   Mean    = hyperparam,
  #   S.D.    = se,
  #   `2.5%`  = hyperparam - qnorm(0.975) * se,
  #   `97.5%` = hyperparam + qnorm(0.975) * se
  # )
  # colnames(param.summ) <- c("Mean", "S.D.", "2.5%", "97.5%")
  # if (length(lambda) > 1) {
  #   rownames(param.summ) <- c("alpha", paste0("lambda[", seq_along(lambda), "]"))
  # } else {
  #   rownames(param.summ) <- c("alpha", "lambda")
  # }
  param.summ <- NULL

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


# w <- x$w
# W <- x$Varw + tcrossprod(w)
# alpha <- get_alpha(x)
# theta <- x$theta
# ystar.minus.alpha <- calc_ystar(x$ipriorKernel, NULL, alpha, theta, w)
#
#
log_q_eta <- function(eta, mod, W, ystar.minus.alpha, w) {
  Hlam <- iprior::.get_Hlam(mod, eta)
  Hlamsq <- iprior::fastSquare(Hlam)
  as.numeric(-0.5 * (
    sum(Hlamsq * W) - 2 * crossprod(ystar.minus.alpha, Hlam) %*% w
  ))
}
#
# eta <- matrix(NA, nrow = 100, ncol = 2)
# eta[1, ] <- theta
# count <- 0
# for (i in seq_len(nrow(eta) - 1)) {
#   eta.star <- rnorm(2, mean = eta[i, ], sd = 0.1)
#   logA <- log_q_eta(eta.star, x$ipriorKernel, W, ystar.minus.alpha, w) -
#     log_q_eta(eta[i, ], x$ipriorKernel, W, ystar.minus.alpha, w)
#   A <- min(exp(logA), 1); print(round(A, 2))
#   if (runif(1) < A) {
#     eta[i + 1, ] <- eta.star
#     count <- count + 1
#   }
#   else {
#     eta[i + 1, ] <- eta[i, ]
#   }
# }

