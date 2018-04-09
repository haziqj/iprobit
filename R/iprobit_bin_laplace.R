#' @export
iprobit_bin_laplace <- function(mod, silent = FALSE, maxit = 100, alpha0 = NULL,
                                theta0 = NULL, w0 = NULL, seed = NULL,
                                stop.crit = 1e-5) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(mod, iprobit.env)
  environment(loop_logical) <- iprobit.env
  maxit <- max(1, maxit)  # cannot have maxit <= 0
  y <- as.numeric(factor(y))  # as.factor then as.numeric to get y = 1, 2, ...

  # Initialise -----------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  if (is.null(alpha0)) alpha0 <- rnorm(1)
  if (is.null(theta0)) theta0 <- rnorm(thetal$n.theta)
  if (is.null(w0)) w0 <- rep(0, n)
  alpha <- alpha0
  theta <- theta0
  w <- w0
  lb <- train.error <- train.brier <- test.error <- test.brier <- rep(NA, maxit)

  # Default optim control list -------------------------------------------------
  control <- list(
    fnscale = -1,
    trace   = ifelse(isTRUE(silent), 0, 1),
    maxit   = max(0, maxit - 1),
    REPORT  = 10,
    factr   = stop.crit / .Machine$double.eps
  )
  # control <- iprior::.update_control(control, control_)

  start.time <- Sys.time()
  res <- optim(c(alpha, theta), lap_bin, object = mod, w = w, trace = TRUE,
               env = environment(), control = control, method = "L-BFGS",
               hessian = TRUE)  # Hlam, w, Varw.inv are written to environment
  lb <- res$value
  Varw <- solve(Varw.inv)
  alpha <- res$par[1]
  theta <- res$par[-1]
  end.time <- Sys.time()
  time.taken <- iprior::as.time(end.time - start.time)

  # Calculate fitted values and error rates ----------------------------------
  eta <- as.numeric(alpha + Hlam %*% w)
  fitted.values <- probs_yhat_error(y, y.levels, eta)
  train.error <- fitted.values$error
  train.brier <- fitted.values$brier
  fitted.test <- NULL
  if (iprior::.is.ipriorKernel_cv(mod)) {
    ystar.test <- calc_ystar(mod, mod$Xl.test, alpha, theta, w)
    fitted.test <- probs_yhat_error(y.test, y.levels, ystar.test)
    test.error[niter + 1] <- fitted.test$error
    test.brier[niter + 1] <- fitted.test$brier
  }

  # Calculate standard errors --------------------------------------------------
  tmp <- iprior::eigenCpp(-res$hessian)
  u <- tmp$val + 1e-9
  V <- tmp$vec
  Fi.inv <- V %*% (t(V) / u)
  se <- sqrt(diag(Fi.inv))
  se[-1] <- iprior::.convert_se(se[-1], theta, mod)  # delta method to convert to parameter s.e.

  # Calculate posterior s.d. and prepare param table ---------------------------
  theta <- matrix(theta, ncol = 2, nrow = length(theta))
  param.full <- theta_to_param.full(theta, alpha, mod)
  param <- param.full[-nrow(param.full), 1]
  param.summ <- cbind(
    param,
    se,
    param - qnorm(0.975) * se,
    param + qnorm(0.975) * se
  )
  colnames(param.summ) <- c("Mean", "S.D.", "2.5%", "97.5%")
  rownames(param.summ) <- get_names(mod, expand = FALSE)

  list(theta = theta, param.full = param.full, param.summ = param.summ, w = w,
       Varw = Varw, lower.bound = as.numeric(na.omit(lb)), niter = res$count[1],
       start.time = start.time, end.time = end.time, time = time.taken,
       fitted.values = fitted.values, test = fitted.test,
       train.error = as.numeric(na.omit(train.error)),
       train.brier = as.numeric(na.omit(train.brier)),
       test.error = as.numeric(na.omit(test.error)),
       test.brier = as.numeric(na.omit(test.brier)), convergence = res$convergence,
       message = res$message)
}

#' @export
lap_bin <- function(mu, object, w0, trace = FALSE, env = NULL) {
  alpha <- mu[1]
  theta <- mu[-1]
  y <- object$y
  y <- as.numeric(factor(y))  # as.factor then as.numeric to get y = 1, 2, ...

  if (missing(w0)) w <- rnorm(object$n)
  else w <- w0

  Hlam <- get_Hlam(object, theta)

  lap.optim <- optim(w, Q_bin, alpha = alpha, Hlam = Hlam, y = y,
                     hessian = TRUE, method = "CG",
                     control = list(fnscale = -1))

  w.tilde <- lap.optim$par
  logdetA <- as.numeric(determinant(-lap.optim$hess)$modulus)

  if (isTRUE(trace)) {
    assign("w", w.tilde, envir = env)
    assign("Varw.inv", -lap.optim$hess, envir = env)
    assign("Hlam", Hlam, envir = env)
  }

  res <- lap.optim$value - logdetA / 2
  res
}

#' @export
Q_bin <- function(w, alpha, Hlam, y) {
  # The Q function for the Laplace method for binary models.
  #
  # Args: w (I-prior random effects), the intercept alpha, the kernel parameters
  # theta, and y (of the form 1, 2, 3, ...)
  #
  # Returns: The Q(w) value.
  f <- alpha + Hlam %*% w
  thing1 <- pnorm(f[y == 2], log.p = TRUE)
  thing0 <- pnorm(-f[y == 1], log.p = TRUE)
  sum(thing1) + sum(thing0) - sum(w ^ 2) / 2
}

# dat <- gen_mixture(100)
# mod <- iprior::kernL(y ~ ., dat, one.lam = TRUE, est.psi = FALSE)
# iprobit_bin_laplace(mod)

