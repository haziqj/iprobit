iprobit_bin_metr <- function(mod, maxit = 20, stop.crit = 1e-5, silent = FALSE,
                             alpha0 = NULL, theta0 = NULL, w0 = NULL,
                             n.samp = 1000, sd.samp = 0.1, thin.samp = 10) {
  # Declare all variables and functions to be used into environment ------------
  iprobit.env <- environment()
  list2env(mod, iprobit.env)
  environment(loop_logical) <- iprobit.env
  maxit <- max(1, maxit)  # cannot have maxit <= 0
  y <- as.numeric(factor(y))  # as.factor then as.numeric to get y = 1, 2, ...
  thin.seq <- seq_len(n.samp)[seq_len(n.samp) %% thin.samp == 0]

  # Initialise -----------------------------------------------------------------
  if (is.null(alpha0)) alpha0 <- rnorm(1)
  if (is.null(theta0)) theta0 <- rnorm(thetal$n.theta)
  if (is.null(w0)) w0 <- rep(0, n)
  alpha <- alpha0
  theta <- theta0
  Hlaml <- Hlamsql <- list(NULL)
  Hlaml[[1]] <- Hlam <- iprior::.get_Hlam(mod, theta)
  Hlamsql[[1]] <- Hlamsq <- iprior::fastSquare(Hlam)
  theta.samp <- matrix(NA, nrow = n.samp, ncol = thetal$n.theta)
  theta.samp[1, ] <- theta
  w <- w0
  Varw <- NA
  log.q <- rep(NA, n.samp)
  log.q[1] <- -0.5 * sum(diag(Hlamsq))
  eta <- as.numeric(alpha + Hlam %*% w)
  niter <- 0
  acc.rate <- lb <- train.error <- train.brier <- test.error <- test.brier <- rep(NA, maxit)
  lb.const <- (n + p - log(n) + p * log(2 * pi)) / 2

  # The variational EM loop ----------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (loop_logical()) {  # see loop_logical() function in iprobit_helper.R
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
    count <- 0
    for (i in seq_len(n.samp - 1)) {
      theta.star <- rnorm(thetal$n.theta, mean = theta.samp[i, ], sd = sd.samp)
      Hlam.star <- iprior::.get_Hlam(mod, theta.star)
      Hlamsq.star <- iprior::fastSquare(Hlam.star)
      log.q.star <- as.numeric(-0.5 * (
        sum(Hlamsq.star * W) - 2 * crossprod(ystar - alpha, Hlam.star) %*% w
      ))
      log.prob.acc <- log.q.star - log.q[i]
      prob.acc <- min(exp(log.prob.acc), 1)
      if (runif(1) < prob.acc) {
        theta.samp[i + 1, ] <- theta.star
        log.q[i + 1] <- log.q.star
        Hlaml[[i + 1]] <- Hlam.star
        Hlamsql[[i + 1]] <- Hlamsq.star
        count <- count + 1
      } else {
        theta.samp[i + 1, ] <- theta.samp[i, ]
        log.q[i + 1] <- log.q[i]
        Hlaml[[i + 1]] <- Hlaml[[i]]
        Hlamsql[[i + 1]] <- Hlamsql[[i]]
      }
    }
    acc.rate[niter] <- count / n.samp
    theta <- apply(theta.samp[thin.seq, ], 2, mean)
    Hlam <- Reduce("+", Hlaml[thin.seq]) / length(Hlaml[thin.seq])
    Hlamsq <- Reduce("+", Hlamsql[thin.seq]) / length(Hlamsql[thin.seq])

    # Update alpha -------------------------------------------------------------
    alpha <- mean(ystar - Hlam %*% w)

    # Calculate lower bound ----------------------------------------------------
    lb[niter + 1] <- lb.const +
      sum(pnorm(eta[y == 2], log.p = TRUE)) +
      sum(pnorm(-eta[y == 1], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus) / 2 - mean(log.q)

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
  param.summ <- rbind(
    alpha = c(alpha, 1 / n, alpha - qnorm(0.975) / n, alpha + qnorm(0.975) / n),
    theta.samp_to_param.summ(theta.samp[thin.seq, ], mod)
  )

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
       test.brier = as.numeric(na.omit(test.brier)), convergence = convergence,
       acc.rate = as.numeric(na.omit(acc.rate)),
       theta.samp = theta.samp[thin.seq, ])
}

theta.samp_to_param.summ <- function(theta.samp, object) {
  # Args: theta.samp is the sample from the (approximate) posterior q(theta),
  # and object is an ipriorKernel object.
  theta <- apply(theta.samp, 2, mean)
  x <- as.data.frame(t(theta.samp))
  ginv.theta.samp <- t(sapply(x, theta_to_coef, object = object))
  ginv.theta <- theta_to_coef(theta, object)
  ginv.theta.sd <- apply(ginv.theta.samp, 2, sd)
  ginv.theta.0025 <- apply(ginv.theta.samp, 2, stats::quantile, probs = 0.025)
  ginv.theta.0975 <- apply(ginv.theta.samp, 2, stats::quantile, probs = 0.975)
  param.summ <- data.frame(
    Mean    = ginv.theta,
    S.D.    = ginv.theta.sd,
    `2.5%`  = ginv.theta.0025,
    `97.5%` = ginv.theta.0975
  )
  colnames(param.summ) <- c("Mean", "S.D.", "2.5%", "97.5%")
  param.summ
}


