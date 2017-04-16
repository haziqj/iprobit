#' @export
fitted.ipriorProbit <- function(x, upper.or.lower = NULL) {
  ystar <- x$ystar
  y.hat <- rep(0, length(x$ystar)); y.hat[ystar >= 0] <- 1
  se.ystar <- iprobitSE(y = y.hat, eta = ystar)

  if (!is.null(upper.or.lower)) {
    if (upper.or.lower == "upper") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 2.241403 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] + 0.03133798 * se.ystar[ystar < 0]
    } else if (upper.or.lower == "lower") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] - 0.03133798 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 2.241403 * se.ystar[ystar < 0]
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
      ystar[ystar >= 0] <- ystar[ystar >= 0] + 2.241403 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] + 0.03133798 * se.ystar[ystar < 0]
    } else if (upper.or.lower == "lower") {
      ystar[ystar >= 0] <- ystar[ystar >= 0] - 0.03133798 * se.ystar[ystar >= 0]
      ystar[ystar < 0] <- ystar[ystar < 0] - 2.241403 * se.ystar[ystar < 0]
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
