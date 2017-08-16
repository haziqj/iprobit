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

#' @export
predict.iprobitMod <- function(object, newdata = list(), y.test = NULL,
                               upper.or.lower = NULL, round.digits = 4, ...) {
  list2env(object, environment())
  list2env(ipriorKernel, environment())
  list2env(model, environment())

  if (length(newdata) == 0) {
    return(cat("No new data supplied. Use fitted() instead."))
  } else {
    if (!is.null(object$formula)) {
      # Model has been fitted using formula interface
      if (is.iprobitData(newdata)) newdata <- as.data.frame(newdata)
      # mf <- model.frame(formula = object$formula, data = newdata)
      tt <- object$ipriorKernel$terms
      Terms <- delete.response(tt)
      xstar <- model.frame(Terms, newdata)
      if (any(colnames(newdata) == yname))
        y.test <- model.extract(model.frame(tt, newdata), "response")
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
    Hl <- .hMatList(x, kernel, intr, no.int, Hurst, intr.3plus,
                    rootkern = FALSE, xstar)

    # Pass to appropriate prediction function ----------------------------------
    if (is.iprobitMod_bin(object)) {
      environment(lambdaExpand_bin) <- environment()
      environment(HlamFn) <- environment()
      lambdaExpand_bin(env = environment(), y = NULL)
      HlamFn(env = environment())
      ystar.new <- as.vector(alpha + (Hlam.mat %*% w))
      names(ystar.new) <- xrownames
      res <- predict_iprobit_bin(y.test, y.levels, ystar.new)
    }
    if (is.iprobitMod_mult(object)) {
      environment(lambdaExpand_mult) <- environment()
      environment(HlamFn_mult) <- environment()
      lambdaExpand_mult(env = environment(), y = NULL)
      HlamFn_mult(env = environment())
      ystar.new <- rep(alpha, each = nrow(Hlam.mat[[1]])) +
        mapply("%*%", Hlam.mat, split(w, col(w)))
      names(ystar.new) <- xrownames
      res <- predict_iprobit_mult(y.test, y.levels, ystar.new)
    }
  }

  res$test.error <- res$train.error
  res$train.error <- NULL
  res
}

brier_score <- function(y, y.hat, prob) {
  prob <- apply(prob, 1, max)
  outcome <- (y == y.hat)
  res <- mean((prob - outcome) ^ 2)
  res
}

predict_iprobit_bin <- function(y, y.levels, ystar) {
  y.hat <- rep(1, length(ystar))
  y.hat[ystar >= 0] <- 2
  y.hat <- factor(y.hat, levels = 1:2)
  levels(y.hat) <- y.levels

  # if (!is.null(upper.or.lower)) {
  #   se.ystar <- iprobitSE(y = y.hat, eta = ystar)
  #   if (upper.or.lower == "upper") {
  #     ystar[ystar >= 0] <- ystar[ystar >= 0] + 2.241403 * se.ystar[ystar >= 0]
  #     ystar[ystar < 0] <- ystar[ystar < 0] - 0.03133798 * se.ystar[ystar < 0]
  #   } else if (upper.or.lower == "lower") {
  #     ystar[ystar >= 0] <- ystar[ystar >= 0] + 0.03133798 * se.ystar[ystar >= 0]
  #     ystar[ystar < 0] <- ystar[ystar < 0] - 2.241403 * se.ystar[ystar < 0]
  #   }
  #   y.hat[ystar >= 0] <- 1
  # }

  p.hat <- pnorm(ystar)
  p.hat <- data.frame(1 - p.hat, p.hat)
  colnames(p.hat) <- y.levels

  error.rate <- mean(as.numeric(y.hat) != as.numeric(y)) * 100
  brier.score <- brier_score(as.numeric(y), as.numeric(y.hat), as.data.frame(p.hat))

  structure(list(y = y.hat, prob = as.data.frame(p.hat),
                 train.error = error.rate, brier.score = brier.score),
            class = "iprobit_predict")
}

predict_iprobit_mult <- function(y, y.levels, ystar) {
  m <- length(y.levels)
  y.hat <- factor(
    apply(ystar, 1, function(x) which(x == max(x))),
    levels = seq_len(m)
  )
  levels(y.hat) <- y.levels

  p.hat <- ystar
  for (i in seq_len(nrow(ystar))) {
    for (j in seq_along(y.levels)) {
      p.hat[i, j] <- EprodPhiZ(ystar[i, j] - ystar[i, seq_along(y.levels)[-j]])
    }
    # p.hat[i, m] <- 1 - sum(p.hat[i, 1:(m - 1)])
  }
  p.hat <- p.hat / matrix(rep(apply(p.hat, 1, sum), m), ncol = m)  # normalise
  colnames(p.hat) <- y.levels

  error.rate <- mean(as.numeric(y.hat) != as.numeric(y)) * 100
  brier.score <- brier_score(as.numeric(y), as.numeric(y.hat), as.data.frame(p.hat))

  structure(list(y = y.hat, prob = as.data.frame(p.hat),
                 train.error = error.rate, brier.score = brier.score),
            class = "iprobit_predict")
}

#' @export
print.iprobitPredict <- function(x, ...) {
  if (!is.null(x$test.error.rate))
    cat("Test error rate:", x$test.error.rate, "%\n")
  else
    cat("Test data not provided.\n")

  cat("\nPredicted classes:\n")
  print(x$y)

  cat("\nPredicted probabilities:\n")
  print(head(x$p, 10), justify = "right")
  if (nrow(x$p > 10)) cat("...")
}

#' @export
print.iprobit_predict <- function(x, rows = 10, dp = 3, ...) {
  if (!is.null(x$train.error)) {
    cat("Training error rate:", decimal_place(x$train.error, dp), "%\n")
    cat("Brier score:", decimal_place(x$brier.score, dp), "\n")
  } else if (!is.nan(x$test.error)) {
    cat("Test error rate:", decimal_place(x$test.error, dp), "%\n")
    cat("Brier score:", decimal_place(x$brier.score, dp), "\n")
  } else {
    cat("Test data not provided.\n")
  }

  cat("\nPredicted classes:\n")
  print(x$y)

  cat("\nPredicted probabilities:\n")
  rows <- min(nrow(x$p), rows)
  tab <- decimal_place(x$p[seq_len(rows), ], dp)
  print(tab)
  if (nrow(x$p) > rows) cat("# ... with", nrow(x$p) - rows, "more rows")
}

# Note: Quantiles for truncated normal distribution are
# qtruncnorm(0.025, a = 0)  # upper tail truncated at zero
# qtruncnorm(0.975, a = 0)
# ## 0.03133798
# ## 2.241403
# lower tail truncated at zero are symmetric opposites of the above.
