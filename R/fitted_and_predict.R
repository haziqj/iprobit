#' @export
fitted.iprobitMod_bin <- function(x, upper.or.lower = NULL, round.digits = 4) {
  ystar <- x$ystar
  y.hat <- rep(0, length(x$ystar)); y.hat[ystar >= 0] <- 1
  se.ystar <- iprobitSE(y = y.hat, eta = ystar)

  if (!is.null(upper.or.lower)) {
    if (upper.or.lower == "upper") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 2.241403 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 0.03133798 * se.ystar[ystar < 0]
    } else if (upper.or.lower == "lower") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 0.03133798 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 2.241403 * se.ystar[ystar < 0]
    }
    y.hat[ystar >= 0] <- 1
  }
  p.hat <- pnorm(ystar)
  p.hat <- cbind(1 - p.hat, p.hat)
  p.hat <- round(p.hat, round.digits)
  colnames(p.hat) <- x$y.lev
  y.hat <- as.factor(y.hat); levels(y.hat) <- x$y.levels

  list(y = y.hat, prob = as.data.frame(p.hat))
}

#' @export
fitted.iprobitMod_mult <- function(x, round.digits = 4) {
  list2env(x, envir = environment())
  n <- length(y)
  y.lev <- levels(y)
  m <- length(y.lev)

  y.hat <- factor(apply(ystar, 1, function(x) which(x == max(x))))
  levels(y.hat) <- y.lev

  probs <- ystar
  for (i in 1:n) {
    for (j in 1:m) {
      probs[i, j] <- EprodPhiZ(ystar[i, j] - ystar[i, (1:m)[-j]])
    }
    # probs[i, m] <- 1 - sum(probs[i, 1:(m - 1)])
  }
  probs <- round(probs, round.digits)
  probs <- probs / matrix(rep(apply(probs, 1, sum), m), ncol = m)  # normalise
  probs <- round(probs, round.digits)
  colnames(probs) <- y.lev

  list(y = y.hat, prob = as.data.frame(probs))
}

#' @export
predict.iprobitMod <- function(object, newdata = list(), y.test = NULL,
                               upper.or.lower = NULL) {
  list2env(object$ipriorKernel, environment())
  list2env(model, environment())
  environment(.lambdaExpand) <- environment()

  if (length(newdata) == 0) {
    return("Use fitted instead")
  } else {
    if (!is.null(object$formula)) {
      # Model has been fitted using formula interface
      if (is.iprobitData(newdata)) newdata <- as.data.frame(newdata)
      mf <- model.frame(formula = object$formula, data = newdata)
      tt <- terms(mf)
      Terms <- delete.response(tt)
      xstar <- model.frame(Terms, newdata)
      y.test <- newdata[, attr(tt, "response")]
      xrownames <- rownames(xstar)
      if (one.lam) {
        xstar <- list(as.matrix(xstar))
      }
    } else {
      if (any(sapply(newdata, is.vector))) {
        newdata <- lapply(newdata, function(x) t(as.matrix(x)))
      }
      xstar <- newdata
      xrownames <- rownames(do.call(cbind, newdata))
    }

    # Define new kernel matrix -------------------------------------------------
    Hl <- .hMatList(x, kernel, intr, no.int, model$Hurst, intr.3plus,
                    rootkern = FALSE, xstar)
    .lambdaExpand(object$lambda, env = environment())
    # if (rootkern) {
    #   Hlam.mat <- object$psi *
    #     Reduce("+", mapply("*", Hl, lambda ^ 2, SIMPLIFY = FALSE))
    #   w.hat <- varyinv(object) %*% (Y - object$alpha)
    # } else {
    Hlam.mat <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
    w.hat <- object$w
    # }

    # Calculate fitted ystar values --------------------------------------------
    ystar <- as.vector(object$alpha + (Hlam.mat %*% w.hat))
    names(ystar) <- xrownames
  }

  if (is.iprobitMod_bin(object)) yp <- predict_iprobit_bin(ystar, y.levels)
  if (is.iprobitMod_mult(object)) yp <- predict_iprobit_mult(ystar, y.levels)

  test.error.rate <- NULL
  if (!is.null(y.test)) {
    test.error.rate <- round(mean(yp$y != y.test) * 100, 2)
  }

  structure(list(y = yp$y, prob = yp$p, test.error.rate = test.error.rate),
            class = "iprobitPredict")
}

predict_iprobit_bin <- function(ystar, y.levels) {
  y.hat <- rep(1, length(ystar))
  y.hat[ystar >= 0] <- 2
  y.hat <- as.factor(y.hat)

  # se.ystar <- iprobitSE(y = y.hat, eta = ystar)
  # if (!is.null(upper.or.lower)) {
  #   if (upper.or.lower == "upper") {
  #     ystar[ystar >= 0] <- ystar[ystar >= 0] + 2.241403 * se.ystar[ystar >= 0]
  #     ystar[ystar < 0] <- ystar[ystar < 0] - 0.03133798 * se.ystar[ystar < 0]
  #   } else if (upper.or.lower == "lower") {
  #     ystar[ystar >= 0] <- ystar[ystar >= 0] + 0.03133798 * se.ystar[ystar >= 0]
  #     ystar[ystar < 0] <- ystar[ystar < 0] - 2.241403 * se.ystar[ystar < 0]
  #   }
  #   y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
  # }

  p.hat <- pnorm(ystar)
  p.hat <- data.frame(1 - p.hat, p.hat)
  y.levels

  colnames(p.hat) <- levels(y.hat) <- y.levels

  list(y.hat = y.hat, p.hat = as.data.frame(p.hat))
}

#' #' @export
#' fitted2 <- function(x, upper.or.lower = NULL) {
#'   ystar <- x$ystar
#'   y.hat <- rep(0, length(x$ystar)); y.hat[ystar >= 0] <- 1
#'   se.ystar <- iprobitSE(y = y.hat, eta = ystar)
#'
#'   if (!is.null(upper.or.lower)) {
#'     if (upper.or.lower == "upper") {
#'       ystar[ystar >= 0] <- ystar[ystar >= 0] + 1.96 * se.ystar[ystar >= 0]
#'       ystar[ystar < 0] <- ystar[ystar < 0] + 1.96 * se.ystar[ystar < 0]
#'     } else if (upper.or.lower == "lower") {
#'       ystar[ystar >= 0] <- ystar[ystar >= 0] - 1.96 * se.ystar[ystar >= 0]
#'       ystar[ystar < 0] <- ystar[ystar < 0] - 1.96 * se.ystar[ystar < 0]
#'     }
#'     y.hat[ystar >= 0] <- 1
#'   }
#'   p.hat <- pnorm(ystar)
#'   y.hat <- as.factor(y.hat); levels(y.hat) <- x$y.levels
#'
#'   list(y = y.hat, prob = p.hat)
#' }



#' #' @export
#' predict2 <- function(object, newdata, upper.or.lower = NULL) {
#'   w <- object$w
#'   lambda <- object$lambda
#'   alpha <- object$alpha
#'
#'   H.tilde <- ikernL(Xl = list(object$X), newdata = list(newdata),
#'                     kernel = object$kernel)[[1]]
#'   class(H.tilde) <- NULL
#'   ystar <- as.numeric(alpha + lambda * H.tilde %*% w)
#'   y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
#'   se.ystar <- iprobitSE(y = y.hat, eta = ystar)
#'
#'   if (!is.null(upper.or.lower)) {
#'     if (upper.or.lower == "upper") {
#'       ystar[ystar >= 0] <- ystar[ystar >= 0] + 1.96 * se.ystar[ystar >= 0]
#'       ystar[ystar < 0] <- ystar[ystar < 0] + 1.96 * se.ystar[ystar < 0]
#'     } else if (upper.or.lower == "lower") {
#'       ystar[ystar >= 0] <- ystar[ystar >= 0] - 1.96 * se.ystar[ystar >= 0]
#'       ystar[ystar < 0] <- ystar[ystar < 0] - 1.96 * se.ystar[ystar < 0]
#'     }
#'     y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
#'   }
#'   p.hat <- pnorm(ystar)
#'   y.hat <- as.factor(y.hat); levels(y.hat) <- object$y.levels
#'
#'   list(y = y.hat, prob = p.hat)
#' }

# Note: Quantiles for truncated normal distribution are
# qtruncnorm(0.025, a = 0)  # upper tail truncated at zero
# qtruncnorm(0.975, a = 0)
# ## 0.03133798
# ## 2.241403
# lower tail truncated at zero are symmetric opposites of the above.



#' @export
predict.iprobitMult <- function(object, X.test, y.test = NULL,
                                upper.or.lower = NULL) {
  newdata <- X.test
  w <- object$w
  n <- nrow(w); m <- ncol(w); n.new <- nrow(newdata)
  lambda <- 0; lambda[1:m] <- object$lambda
  alpha <- 0; alpha[1:m] <- object$alpha
  if (object$kernel == "FBM")
    H.tilde <- iprior::fnH3(x = object$X, y = newdata)
  else
    H.tilde <- iprior::fnH2(x = object$X, y = newdata)
  class(H.tilde) <- NULL
  ystar <- rep(alpha, each = n.new) + rep(lambda, each = n.new) * (H.tilde %*% w)

  y.hat <- apply(ystar, 1, function(x) which(x == max(x)))

  # se.ystar <- iprobitSE(y = y.hat, eta = ystar)
  # if (!is.null(upper.or.lower)) {
  #   if (upper.or.lower == "upper") {
  #     ystar[ystar >= 0] <- ystar[ystar >= 0] + 2.241403 * se.ystar[ystar >= 0]
  #     ystar[ystar < 0] <- ystar[ystar < 0] - 0.03133798 * se.ystar[ystar < 0]
  #   } else if (upper.or.lower == "lower") {
  #     ystar[ystar >= 0] <- ystar[ystar >= 0] + 0.03133798 * se.ystar[ystar >= 0]
  #     ystar[ystar < 0] <- ystar[ystar < 0] - 2.241403 * se.ystar[ystar < 0]
  #   }
  #   y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
  # }
  p.hat <- ystar
  for (i in 1:n.new) {
    for (j in 1:m) {
      p.hat[i, j] <- EprodPhiZ(ystar[i, j] - ystar[i, (1:m)[-j]])
    }
    # probs[i, m] <- 1 - sum(probs[i, 1:(m - 1)])
  }
  p.hat <- round(p.hat, 3)
  p.hat <- p.hat / matrix(rep(apply(p.hat, 1, sum), m), ncol = m)  # normalise
  p.hat <- round(p.hat, 3)
  y.hat <- as.factor(y.hat)
  colnames(p.hat) <- levels(y.hat) <- levels(object$y)

  test.error.rate <- NULL
  if (!is.null(y.test)) {
    test.error.rate <- round(mean(y.hat != y.test) * 100, 2)
  }

  structure(list(y = y.hat, prob = as.data.frame(p.hat),
                 test.error.rate = test.error.rate),
            class = "iprobitMultPredict")
}

#' @export
print.iprobitPredict <- function(x) {
  if (!is.null(x$test.error.rate))
    cat("Test error rate:", x$test.error.rate, "%\n\n")
  else
    cat("Test data not provided.\n\n")
}
