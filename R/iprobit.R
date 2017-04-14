#' @export
iprobitSE <- function(y, eta, thing1 = NULL, thing0 = NULL) {
  if (is.null(thing1) | is.null(thing0)) {
    thing1 <- exp(  # phi(eta) / Phi(eta)
      dnorm(eta[y == 1], log = TRUE) - pnorm(eta[y == 1], log.p = TRUE)
    )
    thing0 <- -exp(  # -1 * {phi(eta) / Phi(-eta)}
      dnorm(eta[y == 0], log = TRUE) - pnorm(-eta[y == 0], log.p = TRUE)
    )
  }

  # Posterior variance of ystar ------------------------------------------------
  var.ystar <- rep(NA, length(y))
  # 1 - eta * phi(eta) / Phi(eta) - (phi(eta) / Phi(eta)) ^ 2
  var.ystar[y == 1] <- 1 - eta[y == 1] * thing1 + (thing1 ^ 2)
  # 1 - eta * (-1) * {phi(eta) / Phi(-eta)} - (phi(eta) / Phi(-eta)) ^ 2
  var.ystar[y == 0] <- 1 - eta[y == 0] * thing0 + (thing0 ^ 2)
  sqrt(var.ystar)
}

#' @export
iprobit <- function(y, ..., kernel = "Canonical", maxit = 1000, stop.crit = 1e-5,
                    silent = FALSE, interactions = NULL, alpha0 = rnorm(1),
                    lambda0 = abs(rnorm(1)), w0 = rep(0, n)) {
  y.tmp <- checkLevels(y)
  y <- y.tmp$y
  y.levels <- y.tmp$levels

  # Prepare kernel matrices ----------------------------------------------------
  Xl <- list(...)
  iprobit.kernel <- ikernL(Xl, NULL, kernel, NULL)
  H <- iprobit.kernel[[1]]  # CHANGE THIS FOR FUTURE UPDATES. NOW ONLY SUPPORT
                            # SINGLE LAMBDA

  H2 <- H %*% H
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
  n <- length(y)

  # Set up parameter results ---------------------------------------------------
  lower.bound <- lambda <- rep(NA, maxit)
  error.rate <- alpha <- rep(0, maxit)
  w <- matrix(NA, ncol = n, nrow = maxit)
  ystar <- matrix(NA, ncol = n, nrow = maxit)

  # Initialise -----------------------------------------------------------------
  lambda[1] <- lambda0
  lambda2 <- lambda0 ^ 2
  alpha[1] <- alpha0
  w[1, ] <- w0
  niter <- 1
  lb.const <- (n + 2 - log(n)) / 2 + log(2 * pi)

  start.time <- Sys.time()
  for (t in 1:(maxit - 1)) {
    # Update ystar -------------------------------------------------------------
    eta <- as.numeric(alpha[t] + lambda[t] * H %*% w[t, ])
    thing <- rep(NA, n)
    thing1 <- exp(  # phi(eta) / Phi(eta)
      dnorm(eta[y == 1], log = TRUE) - pnorm(eta[y == 1], log.p = TRUE)
    )
    thing0 <- -exp(  # -1 * {phi(eta) / Phi(-eta)}
      dnorm(eta[y == 0], log = TRUE) - pnorm(-eta[y == 0], log.p = TRUE)
    )
    thing[y == 1] <- thing1
    thing[y == 0] <- thing0
    ystar[t + 1, ] <- eta + thing

    # Update w -----------------------------------------------------------------
    A <- lambda2 * H2 + diag(1, n)
    a <- as.numeric(lambda[t] * crossprod(H, ystar[t + 1, ] - alpha[t]))
    eigenA <- eigen(A)
    V <- eigenA$vec
    u <- eigenA$val + 1e-8  # ensure positive eigenvalues
    uinv.Vt <- t(V) / u
    w[t + 1, ] <- as.numeric(crossprod(a, V) %*% uinv.Vt)
    Varw <- V %*% uinv.Vt
    W <- Varw + tcrossprod(w[t + 1, ])

    # Update lambda ------------------------------------------------------------
    ct <- sum(H2 * W)
    d <- as.numeric(crossprod(ystar[t + 1, ] - alpha[t], H) %*% w[t + 1, ])
    lambda[t + 1] <- d / ct
    lambda2 <- 1 / ct + (d / ct) ^ 2

    # Update alpha -------------------------------------------------------------
    alpha[t + 1] <- mean(ystar[t + 1, ] - lambda[t + 1] * H %*% w[t + 1, ])

    # Lower bound --------------------------------------------------------------
    lower.bound[t + 1] <- lb.const +
      sum(pnorm(eta[y == 1], log.p = TRUE)) + sum(pnorm(-eta[y == 0], log.p = TRUE)) -
      (sum(diag(W)) + determinant(A)$modulus + log(ct)) / 2

    # Running fit --------------------------------------------------------------
    tmp <- rep(0, n)
    tmp[ystar[t + 1, ] >= 0] <- 1
    error.rate[t] <- sum(tmp != y) / n * 100

    lb.diff <- abs(lower.bound[t + 1] - lower.bound[t])
    if (!is.na(lb.diff) && (lb.diff < stop.crit)) break
    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, t)
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  # Calculate standard errors from posterior variance --------------------------
  se.lambda <- sqrt(1 / ct)
  se.alpha <- sqrt(1 / n)
  se.ystar <- iprobitSE(y = y, eta = eta, thing1 = thing1, thing0 = thing0)

  # Close function -------------------------------------------------------------

  if (!silent) {
    close(pb)
    if (niter == maxit) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  res <- list(ystar = ystar[niter, ], w = w[niter, ], lambda = lambda[niter],
              alpha = alpha[niter], lower.bound = lower.bound, kernel = kernel,
              X = Xl[[1]], y = y, error.rate = error.rate,
              se = c(se.alpha, se.lambda),
              se.ystar = se.ystar, y.levels = y.levels, Varw = Varw,
              start.time = start.time, end.time = end.time, time = time.taken,
              call = match.call(), stop.crit = stop.crit, niter = niter,
              maxit = maxit)
  class(res) <- "ipriorProbit"
  res
}

#' @export
ipriorProbitPrintAndSummary <- function(x) {
  y.hat <- fitted.ipriorProbit(x)$y
  y <- as.factor(x$y); levels(y) <- x$y.levels
  train.error.rate <- format(round(mean(y.hat != y) * 100, 2))

  # Calculate 95% credibility interval for error rate --------------------------
  y.hat.upper <- fitted(x, "upper")$y
  train.error.rate.upper <- format(round(mean(y.hat.upper != x$y) * 100, 2))
  y.hat.lower <- fitted(x, "lower")$y
  train.error.rate.lower <- format(round(mean(y.hat.lower != x$y) * 100, 2))

  list(
    train.error.rate = train.error.rate,
    error.band = c(train.error.rate.lower, train.error.rate.upper),
    lb = as.numeric(logLik(x))
  )
}

#' @export
print.ipriorProbit <- function(x, newdata = NULL, testdata = NULL) {
  tmp <- ipriorProbitPrintAndSummary(x)
  train.error.rate <- tmp$train.error.rate
  cat("Training error rate:", train.error.rate, "%")
  if (!is.null(newdata) && !is.null(testdata)) {
    y.hat <- predict(x, newdata)$y
    test.error.rate <- format(round(mean(y.hat != testdata) * 100, 2))
    cat("\nTest error rate:", test.error.rate, "%")
  }
}


#' @export
summary.ipriorProbit <- function(x) {
  tmp <- ipriorProbitPrintAndSummary(x)
  train.error.rate <- tmp$train.error.rate

  post.mean <- c(x$alpha, x$lambda)
  se <- x$se  # only for lambda alpha and lambda
  tab <- cbind(
    Mean    = round(post.mean, digits = 4),
    S.E.    = round(se, digits = 4),
    "2.5%"  = round(post.mean - 1.96 * se, digits = 4),
    "97.5%" = round(post.mean + 1.96 * se, digits = 4)
  )
  rownames(tab) <- c("alpha", paste0("lambda", 1:length(x$lambda)))

  res <- list(tab = tab, call = x$call, kernel = x$kernel, maxit = x$maxit,
              stop.crit = x$stop.crit, niter = x$niter, lb = tmp$lb,
              error = train.error.rate, error.band = x$error.band)
  class(res) <- "iprobitSummary"
  res
}

#' @export
print.iprobitSummary <- function(x) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nRKHS used:", x$kernel, "\n\n")
  print(x$tab)
  if (x$niter == x$maxit) {
    cat("\nConvergence criterion not met. ")
  } else {
    cat("\nConverged to within", x$stop.crit, "tolerance. ")
  }
  cat("No. of iterations:", x$niter)
  # cat("\nModel classification error rate (%):", x$error, "within [",
  #     x$error.band[1], ",", x$error.band[2], "]")
  cat("\nModel classification error rate (%):", x$error)
  cat("\nVariational lower bound:", x$lb)
  cat("\n\n")
}
