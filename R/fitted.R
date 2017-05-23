#' @export
fitted.ipriorProbit <- function(x, upper.or.lower = NULL) {
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
  y.hat <- as.factor(y.hat); levels(y.hat) <- x$y.levels

  list(y = y.hat, prob = p.hat)
}

#' @export
fitted2 <- function(x, upper.or.lower = NULL) {
  ystar <- x$ystar
  y.hat <- rep(0, length(x$ystar)); y.hat[ystar >= 0] <- 1
  se.ystar <- iprobitSE(y = y.hat, eta = ystar)

  if (!is.null(upper.or.lower)) {
    if (upper.or.lower == "upper") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 1.96 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] + 1.96 * se.ystar[ystar < 0]
    } else if (upper.or.lower == "lower") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] - 1.96 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 1.96 * se.ystar[ystar < 0]
    }
    y.hat[ystar >= 0] <- 1
  }
  p.hat <- pnorm(ystar)
  y.hat <- as.factor(y.hat); levels(y.hat) <- x$y.levels

  list(y = y.hat, prob = p.hat)
}

#' @export
predict.ipriorProbit <- function(object, newdata, upper.or.lower = NULL) {
  w <- object$w
  lambda <- object$lambda
  alpha <- object$alpha

  H.tilde <- ikernL(Xl = list(object$X), newdata = list(newdata),
                    kernel = object$kernel)[[1]]
  class(H.tilde) <- NULL
  ystar <- as.numeric(alpha + lambda * H.tilde %*% w)
  y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
  se.ystar <- iprobitSE(y = y.hat, eta = ystar)

  if (!is.null(upper.or.lower)) {
    if (upper.or.lower == "upper") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 2.241403 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 0.03133798 * se.ystar[ystar < 0]
    } else if (upper.or.lower == "lower") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 0.03133798 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 2.241403 * se.ystar[ystar < 0]
    }
    y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
  }
  p.hat <- pnorm(ystar)
  y.hat <- as.factor(y.hat); levels(y.hat) <- object$y.levels

  list(y = y.hat, prob = p.hat)
}

#' @export
predict2 <- function(object, newdata, upper.or.lower = NULL) {
  w <- object$w
  lambda <- object$lambda
  alpha <- object$alpha

  H.tilde <- ikernL(Xl = list(object$X), newdata = list(newdata),
                    kernel = object$kernel)[[1]]
  class(H.tilde) <- NULL
  ystar <- as.numeric(alpha + lambda * H.tilde %*% w)
  y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
  se.ystar <- iprobitSE(y = y.hat, eta = ystar)

  if (!is.null(upper.or.lower)) {
    if (upper.or.lower == "upper") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 1.96 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] + 1.96 * se.ystar[ystar < 0]
    } else if (upper.or.lower == "lower") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] - 1.96 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 1.96 * se.ystar[ystar < 0]
    }
    y.hat <- rep(0, nrow(newdata)); y.hat[ystar >= 0] <- 1
  }
  p.hat <- pnorm(ystar)
  y.hat <- as.factor(y.hat); levels(y.hat) <- object$y.levels

  list(y = y.hat, prob = p.hat)
}

# Note: Quantiles for truncated normal distribution are
# qtruncnorm(0.025, a = 0)  # upper tail truncated at zero
# qtruncnorm(0.975, a = 0)
# ## 0.03133798
# ## 2.241403
# lower tail truncated at zero are symmetric opposites of the above.

#' @export
fitted.iprobitMult <- function(x, round.digits = 3) {
  list2env(x, envir = environment())
  n <- length(y)
  y.lev <- levels(y)
  m <- length(y.lev)
  nm <- n * m
  p <- ncol(X)

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
  colnames(probs) <- y.lev

  list(y.hat = y.hat, probs = as.data.frame(probs))
}

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
print.iprobitMultPredict <- function(x) {
  cat("Test error rate:", x$test.error.rate, "%")
}