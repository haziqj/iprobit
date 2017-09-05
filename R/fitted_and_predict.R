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
fitted.iprobitMod <- function(object, quantiles = TRUE, n.samp = 100,
                              transform = function(x) x, ...) {
  if (isTRUE(quantiles)) {
    if (is.iprobitMod_bin(object)) {
      return(predict_quant(object, n.samp, transform))
    }
    if (is.iprobitMod_mult(object)) {
      return(1)
    }
  } else {
    return(object$fitted.values)
  }

}

predict_quant <- function(object, n.samp, transform, Hl = NULL, y = NULL) {
  tmp <- sample_prob_bin(object, n.samp, Hl, y)
  phats <- convert_prob(tmp$phat.samp, transform)
  names(phats) <- object$ipriorKernel$y.levels
  res <- list(
    prob        = quantile_prob(phats),
    train.error = stats::quantile(tmp$error.samp, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)),
    brier.score = stats::quantile(tmp$brier.samp, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  )
  class(res) <- "iprobit_predict_quant"
  res
}

#' @export
print.iprobit_predict_quant <- function(x, rows = 5, dp = 3, ...) {
  rows.act <- nrow(x$prob[[1]])
  rows <- min(rows.act, rows)
  y.levels <- names(x$prob)

  if (!is.null(x$train.error)) {
    error <- decimal_place(x$train.error, dp)
    brier <- decimal_place(x$brier.score, dp)
  } else if (!is.nan(x$test.error)) {
    # cat("Test error rate:", decimal_place(x$test.error, dp), "%\n")
    # cat("Brier score:", decimal_place(x$brier.score, dp), "\n")
  } else {
    # cat("Test data not provided.\n")
  }

  tab <- rbind("Training error (%)" = error, "Brier score" = brier)
  print(as.data.frame(tab))

  for (j in seq_along(x$prob)) {
    cat("\nPredicted probabilities for Class =", y.levels[j], "\n")
    tab <- decimal_place(x$prob[[j]][seq_len(rows), ], dp)
    print(tab)
    if (rows.act > rows) cat("# ... with", rows.act - rows, "more rows\n")
  }
}

sample_prob_bin <- function(object, n.samp, Hl.new = NULL, y = NULL) {
  # Initialise -----------------------------------------------------------------
  alpha.samp <- rnorm(n.samp, object$alpha, object$se.alpha)
  lambda.samp <- mvtnorm::rmvnorm(n.samp, object$lambda, diag(object$se.lambda ^ 2))
  w.samp <- mvtnorm::rmvnorm(n.samp, object$w, object$Varw)  # n.samp x n matrix
  phat.samp <- list()
  error.samp <- brier.samp <- rep(NA, n.samp)
  if (is.null(y)) y <- object$ipriorKernel$Y

  # Sampling -------------------------------------------------------------------
  for (i in 1:n.samp) {
    alpha <- alpha.samp[i]
    lambda <- lambda.samp[i, ]
    w <- w.samp[i, ]
    ystar.samp <- calc_ystar_bin(object, Hl.new = Hl.new, alpha.new = alpha,
                                 lambda.new = lambda, w.new = w)
    tmp <- predict_iprobit_bin(y, object$ipriorKernel$y.levels, ystar.samp)
    phat.samp[[i]] <- tmp$prob
    error.samp[i] <- tmp$train.error
    brier.samp[i] <- tmp$brier.score
  }

  list(phat.samp = phat.samp, error.samp = error.samp, brier.samp = brier.samp)
}

calc_ystar_bin <- function(object, Hl.new = NULL, alpha.new = NULL,
                           lambda.new = NULL, w.new = NULL) {
  list2env(object, environment())
  list2env(ipriorKernel, environment())
  list2env(model, environment())
  environment(lambdaExpand_bin) <- environment()
  environment(HlamFn) <- environment()
  if (!is.null(Hl.new)) Hl <- Hl.new
  if (!is.null(alpha.new)) alpha <- alpha.new
  if (!is.null(lambda.new)) lambda <- lambda.new
  if (!is.null(w.new)) w <- w.new
  lambdaExpand_bin(env = environment(), y = NULL)
  HlamFn(env = environment())
  as.vector(alpha + (Hlam.mat %*% w))
}

convert_prob <- function(phat, transform = function(x) x) {
  # Converts list of length n.samp containing random samples of the
  # probabilities. This function converts it into a list of length no.classes so
  # that each column is now the random sample of the probability of that data
  # point. The transform option transform the probabilities, say to log odds.
  no.classes <- ncol(phat[[1]])
  res <- list()
  for (i in seq_len(no.classes)) {
    res[[i]] <- do.call(cbind, lapply(lapply(phat, `[`, i), transform))
  }
  res
}

quantile_prob <- function(phats) {
  # Get quantiles for lists returned from convert_prob
  for (i in seq_along(phats)) {
    phats[[i]] <- t(apply(phats[[i]], 1, stats::quantile,
                          probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
  }
  lapply(phats, as.data.frame)
}

#' @export
predict.iprobitMod <- function(object, newdata = list(), y.test = NULL,
                               quantiles = TRUE, n.samp = 100,
                               transform = function(x) x, ...) {
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
      if (isTRUE(quantiles)) {
        1 + 1
        return(predict_quant(object, n.samp, transform, Hl, y.test))
      } else {
        ystar.new <- calc_ystar_bin(object, Hl.new = Hl)
        names(ystar.new) <- xrownames
        res <- predict_iprobit_bin(y.test, y.levels, ystar.new)
      }
    }
    if (is.iprobitMod_mult(object)) {
      environment(lambdaExpand_mult) <- environment()
      environment(HlamFn_mult) <- environment()
      lambdaExpand_mult(env = environment(), y = NULL)
      HlamFn_mult(env = environment())
      ystar.new <- rep(alpha, each = nrow(Hlam.mat[[1]])) +
        mapply("%*%", Hlam.mat, split(w, col(w)))
      names(ystar.new) <- xrownames
      res <- predict_iprobit_mult(y.test, y.levels, ystar.new, sd.shift)
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

  p.hat <- pnorm(ystar)
  p.hat <- data.frame(1 - p.hat, p.hat)
  colnames(p.hat) <- y.levels

  error.rate <- mean(as.numeric(y.hat) != as.numeric(y)) * 100
  brier.score <- brier_score(as.numeric(y), as.numeric(y.hat), as.data.frame(p.hat))

  structure(list(y = y.hat, prob = as.data.frame(p.hat),
                 train.error = error.rate, brier.score = brier.score),
            class = "iprobit_predict")
}

predict_iprobit_mult <- function(y, y.levels, ystar, shift = 0) {
  m <- length(y.levels)

  p.hat <- ystar
  for (i in seq_len(nrow(ystar))) {
    for (j in seq_along(y.levels)) {
      p.hat[i, j] <- EprodPhiZ(ystar[i, j] - ystar[i, seq_along(y.levels)[-j]] + shift)
    }
    # p.hat[i, m] <- 1 - sum(p.hat[i, 1:(m - 1)])
  }
  p.hat <- p.hat / matrix(rep(apply(p.hat, 1, sum), m), ncol = m)  # normalise
  colnames(p.hat) <- y.levels

  if (shift == 0) tmp <- ystar
  else tmp <- p.hat
  y.hat <- factor(
    apply(tmp, 1, function(x) which(x == max(x))),
    levels = seq_len(m)
  )
  levels(y.hat) <- y.levels

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
  rows <- min(nrow(x$p), rows)

  cat("\nPredicted classes:\n")
  y.toprint <- x$y[seq_len(rows)]
  y.output <- capture.output(y.toprint)
  y.levels.output <- y.output[length(y.output)]
  y.levels.output <- gsub(" Levels:", "Levels:", y.levels.output)
  y.output <- y.output[-length(y.output)]
  cat(y.output)
  if (nrow(x$p) > rows) cat(" ...")
  cat("\n")
  cat(y.levels.output, "\n")

  cat("\nPredicted probabilities:\n")
  tab <- decimal_place(x$p[seq_len(rows), ], dp)
  print(tab)
  if (nrow(x$p) > rows) cat("# ... with", nrow(x$p) - rows, "more rows")
}
