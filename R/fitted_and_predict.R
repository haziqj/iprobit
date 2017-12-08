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

# For posterity: Both fitted and predict use the same function to output
# predicted values, probabilities, error rate and brier score. This function is
# predict_iprobit_x(). This function is also used in the variational algorithm
# to get $fitted.values.
#
# predict_iprobit_x() arguments are the ystar values, y/y.test and y.levels. In
# the variational algorithm the ystar are calculated sequentially as part of the
# algorithm. In fitted(), ystar are obtained from object$ystar. In predict(),
# ystar are calculated from values of Hl, lambda, alpha and w. There is a helper
# function calc_ystar_x() to calculate these.
#
# To calculate the quantiles of the predicted values, values of lambda, alpha
# and w are simulated from their posterior distribution, and for each sample a
# value of ystar is calculated and predicted values obtained through
# predict_iprobit_x(). Quantiles are calculated using the function
# predict_quant() for both iprobitMod_x type objects. It uses a helper function
# sample_prob_x() to do the sampling.
#
# predict() and fitted() return an iprobit_predict object, while quantiles
# return an iprobit_predict_quant object. This allows different print methods
# for each.

#' @export
fitted.iprobitMod <- function(object, quantiles = FALSE, n.samp = 100,
                              transform = identity, raw = FALSE, ...) {
  if (isTRUE(quantiles)) {
    return(predict_quant(object, n.samp, transform, raw = raw, type = "train"))
  } else {
    res <- object$fitted.values
    res$prob <- transform(res$prob)
    return(res)
  }
}

#' @export
predict.iprobitMod <- function(object, newdata = list(), y.test = NULL,
                               quantiles = FALSE, n.samp = 100,
                               transform = identity, raw = FALSE, ...) {
  if (length(newdata) == 0) {
    return(cat("No new data supplied. Use fitted() instead."))
  }
  if (!is.null(object$ipriorKernel$formula)) {
    if (is.iprobitData(newdata)) newdata <- as.data.frame(newdata)
    tmp <- iprior::.terms_to_xy(object$ipriorKernel, newdata)
    y.test <- tmp$y
    xstar <- tmp$Xl
    xrownames <- rownames(newdata)
  } else {
    if (any(sapply(newdata, is.vector))) {
      newdata <- lapply(newdata, as.matrix)
    }
    xstar <- newdata
    xrownames <- rownames(do.call(cbind, newdata))
  }

  res <- predict_iprobit(object, xstar, y.test, quantiles = quantiles,
                         n.samp = n.samp, transform = transform, raw = raw)
  # names(res$y) <- xrownames
  res
}

brier_score <- function(y, y.hat, prob) {
  prob <- apply(prob, 1, max)
  outcome <- (y == y.hat)
  res <- mean((prob - outcome) ^ 2)
  res
}

predict_iprobit <- function(object, xstar, y.test, quantiles = FALSE,
                            n.samp = 100, transform = identity, raw = FALSE) {
  # Args: object is ipriorMod_x.
  if (!isTRUE(quantiles)) {
    ystar.new <- calc_ystar(object$ipriorKernel, xstar, get_alpha(object),
                            get_theta(object), object$w)
    res <- probs_yhat_error(y.test, object$ipriorKernel$y.levels, ystar.new,
                            type = "test")
  } else {
    res <- predict_quant(object, n.samp, transform, raw, xstar, y.test,
                         type = "test")
  }
  res
}

calc_ystar <- function(object, xstar, alpha, theta, w) {
  # Args: An ipriorKernel object.
  #
  # Returns: For binary models, a vector of ystar. For multinomial models, this
  # is a matrix with columns representing the number of classes.
  if (is.null(xstar)) Hlam.new <- get_Hlam(object, theta)
  else Hlam.new <- get_Htildelam(object, theta, xstar)
  if (is.iprobit_bin(object)) {
    return(as.numeric(alpha + Hlam.new %*% w))
  } else {
    res <- NULL
    m <- get_m(object)
    if (length(alpha) == 1) alpha <- rep(alpha, m)
    return(rep(alpha, each = nrow(Hlam.new[[1]])) +
             mapply("%*%", Hlam.new, split(w, col(w))))
  }
}

probs_yhat_error <- function(y, y.levels, ystar, type = c("train", "test")) {
  m <- length(y.levels)

  if (m == 2) {
    # Binary model -------------------------------------------------------------
    p.hat <- pnorm(ystar)
    p.hat <- data.frame(1 - p.hat, p.hat)

    y.hat <- rep(1, length(ystar))
    y.hat[ystar >= 0] <- 2
    y.hat <- factor(y.hat, levels = 1:2)
  } else {
    # Multinomial model --------------------------------------------------------
    p.hat <- ystar
    for (i in seq_len(nrow(ystar))) {
      for (j in seq_along(y.levels)) {
        p.hat[i, j] <- EprodPhiZ(ystar[i, j] - ystar[i, seq_along(y.levels)[-j]])
      }
      # p.hat[i, m] <- 1 - sum(p.hat[i, 1:(m - 1)])
    }
    p.hat <- p.hat / matrix(rep(apply(p.hat, 1, sum), m), ncol = m)  # normalise

    y.hat <- factor(
      apply(p.hat, 1, function(x) which(x == max(x))),
      levels = seq_len(m)
    )
  }
  colnames(p.hat) <- y.levels
  levels(y.hat) <- y.levels

  error.rate <- mean(as.numeric(y.hat) != as.numeric(y)) * 100
  brier.score <- brier_score(as.numeric(y), as.numeric(y.hat), as.data.frame(p.hat))

  structure(list(y = y.hat, prob = as.data.frame(p.hat),
                 type = match.arg(type, c("train", "test")),
                 error.rate = error.rate, brier.score = brier.score),
            class = "iprobitPredict")
}

#' @export
predict_quant <- function(object, n.samp, transform = identity, raw = FALSE,
                          xstar = NULL, y = NULL, type = c("train", "test")) {
  # Helper function to get the quantiles of the probabilities, error rate and
  # brier score. This is done by sampling from the posterior of the parameters.
  # An option raw = TRUE returns the sampling values instead.
  #
  # Args: An iprobitMod_x object.
  if (is.iprobitMod_bin(object)) {
    tmp <- sample_prob_bin(object, n.samp, xstar, y)
  }
  if (is.iprobitMod_mult(object)) {
    # tmp <- sample_prob_mult(object, n.samp, Hl, y)
  }

  phats <- convert_prob(tmp$phat.samp, transform)
  names(phats) <- object$ipriorKernel$y.levels
  if (isTRUE(raw)) {
    return(list(
      prob = phats, train.error = tmp$error.samp, brier.score = tmp$brier.samp
    ))
  } else {
    return(structure(list(
      prob        = quantile_prob(phats),
      error.rate  = stats::quantile(tmp$error.samp, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      brier.score = stats::quantile(tmp$brier.samp, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      type        = match.arg(type, c("train", "test"))
    ), class = "iprobitPredict_quant"))
  }
}

# calc_ystar <- function(object, xstar, alpha, theta, w)

sample_prob_bin <- function(object, n.samp, xstar = NULL, y = NULL) {
  # Helper function to sample probabilities, error rates and brier scores for
  # binary models.
  #
  # Args: An iprobitMod_x object.
  alpha.samp <- sample_alpha(n.samp, get_alpha(object), get_sd_alpha(object))
  if (grepl("Closed-form", object$est.method)) {
    theta.samp <- sample_lambda(n.samp, get_lambda(object),
                                get_sd_lambda(object))
  } else {
    stop("Not implemented yet.")
  }
  w.samp <- sample_w(n.samp, object$w, object$Varw)  # n.samp x n matrix
  phat.samp <- list()
  error.samp <- brier.samp <- rep(NA, n.samp)
  if (is.null(y)) y <- as.numeric(factor(object$ipriorKernel$y))

  for (i in seq_len(n.samp)) {
    alpha <- alpha.samp[i]
    theta <- hyperparam_to_theta(theta.samp[i, ])
    w <- w.samp[i, ]
    ystar.samp <- calc_ystar(object$ipriorKernel, xstar, alpha, theta, w)
    tmp <- probs_yhat_error(y, object$ipriorKernel$y.levels, ystar.samp)
    phat.samp[[i]] <- tmp$prob
    error.samp[i] <- tmp$error
    brier.samp[i] <- tmp$brier
  }

  list(phat.samp = phat.samp, error.samp = error.samp, brier.samp = brier.samp)
}

sample_alpha <- function(n.samp, mean, sd) {
  rnorm(n.samp, mean, sd)
}

sample_lambda <- function(n.samp, mean, sd) {
  if (length(mean) > 1) {
    return(mvtnorm::rmvnorm(n.samp, mean, diag(sd ^ 2)))
  } else {
    return(matrix(log(rnorm(n.samp, mean, sd))))  # note the log(lambda)
  }
}

sample_w <- function(n.samp, mean, variance) {
  mvtnorm::rmvnorm(n.samp, mean, variance)
}

#' calc_ystar <- function(object, xstar) {
#'   # Helper function to calculate ystar values from Hl, lambda, alpha and w. If
#'   # values are not supplied in the argument, these are obtained from object.
#'   #
#'   # Args: iprobitMod object.
#'   #
#'   # Returns:
#'   list2env(object, environment())
#'   list2env(ipriorKernel, environment())
#'   list2env(model, environment())
#'   if (!is.null(Hl.new)) Hl <- Hl.new
#'   if (!is.null(alpha.new)) alpha <- alpha.new
#'   if (!is.null(lambda.new)) lambda <- lambda.new
#'   if (!is.null(w.new)) w <- w.new
#'
#'   if (is.iprobitMod_bin(object)) {
#'     environment(lambdaExpand_bin) <- environment()
#'     environment(HlamFn) <- environment()
#'     lambdaExpand_bin(env = environment(), y = NULL)
#'     HlamFn(env = environment())
#'     return(as.vector(alpha + (Hlam.mat %*% w)))
#'   }
#'   if (is.iprobitMod_mult(object)) {
#'     # environment(lambdaExpand_mult) <- environment()
#'     # environment(HlamFn_mult) <- environment()
#'     # lambdaExpand_mult(env = environment(), y = NULL)
#'     # HlamFn_mult(env = environment())
#'     # return(rep(alpha, each = nrow(Hlam.mat[[1]])) +
#'     #          mapply("%*%", Hlam.mat, split(w, col(w))))
#'   }
#' }










# predict_iprobit_mult <- function(y, y.levels, ystar, type = c("train", "test")) {
#   m <- length(y.levels)
#
#   p.hat <- ystar
#   for (i in seq_len(nrow(ystar))) {
#     for (j in seq_along(y.levels)) {
#       p.hat[i, j] <- EprodPhiZ(ystar[i, j] - ystar[i, seq_along(y.levels)[-j]])
#     }
#     # p.hat[i, m] <- 1 - sum(p.hat[i, 1:(m - 1)])
#   }
#   p.hat <- p.hat / matrix(rep(apply(p.hat, 1, sum), m), ncol = m)  # normalise
#   colnames(p.hat) <- y.levels
#
#   tmp <- p.hat
#   y.hat <- factor(
#     apply(tmp, 1, function(x) which(x == max(x))),
#     levels = seq_len(m)
#   )
#   levels(y.hat) <- y.levels
#
#   error.rate <- mean(as.numeric(y.hat) != as.numeric(y)) * 100
#   brier.score <- brier_score(as.numeric(y), as.numeric(y.hat), as.data.frame(p.hat))
#
#   type <- match.arg(type, c("train", "test"))
#
#   structure(list(y = y.hat, prob = as.data.frame(p.hat), type = type,
#                  train.error = error.rate, brier.score = brier.score),
#             class = "iprobit_predict")
# }




#' @export
sample_prob_mult <- function(object, n.samp, Hl.new = NULL, y = NULL) {
  # Helper function to sample probabilities, error rates and brier scores for
  # multinomial models.
  alpha.samp <- lambda.samp <- w.samp <- phat.samp <- list()
  m <- object$ipriorKernel$m
  for (j in seq_len(m)) {
    # First sample per class
    alpha.samp[[j]] <- rnorm(n.samp, object$alpha[j], object$se.alpha)
    if (nrow(object$lambda) > 1) {
      lambda.samp[[j]] <- mvtnorm::rmvnorm(n.samp, object$lambda[, j],
                                           diag(object$se.lambda[, j] ^ 2))
    } else {
      lambda.samp[[j]] <- matrix(rnorm(n.samp, object$lambda[, j],
                                       object$se.lambda[, j]))
    }
    w.samp[[j]] <- mvtnorm::rmvnorm(n.samp, object$w[, j], object$Varw[[j]])
  }
  # Then combine each sample into a list of length n.samp (makes it easier in
  # the calculation of the probabilities etc.)
  alpha.samp <- mapply(function(x) unlist(lapply(alpha.samp, `[`, x)),
                       seq_len(n.samp), SIMPLIFY = FALSE)
  lambda.samp <- mapply(function(x) matrix(unlist(lapply(lambda.samp, `[`, x, )), ncol = m),
                        seq_len(n.samp), SIMPLIFY = FALSE)
  w.samp <- mapply(function(x) matrix(unlist(lapply(w.samp, `[`, x, )), ncol = m),
                   seq_len(n.samp), SIMPLIFY = FALSE)
  phat.samp <- list()
  error.samp <- brier.samp <- rep(NA, n.samp)
  if (is.null(y)) y <- object$ipriorKernel$Y

  for (i in seq_len(n.samp)) {
    alpha <- alpha.samp[[i]]
    lambda <- lambda.samp[[i]]
    w <- w.samp[[i]]
    ystar.samp <- calc_ystar(object, Hl.new = Hl.new, alpha.new = alpha,
                             lambda.new = lambda, w.new = w)
    tmp <- predict_iprobit_mult(y, object$ipriorKernel$y.levels, ystar.samp)
    phat.samp[[i]] <- tmp$prob
    error.samp[i] <- tmp$train.error
    brier.samp[i] <- tmp$brier.score
  }

  list(phat.samp = phat.samp, error.samp = error.samp, brier.samp = brier.samp)
}



#' @export
convert_prob <- function(phat, transform = identity) {
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

#' @export
quantile_prob <- function(phats) {
  # Get quantiles for lists returned from convert_prob()
  for (i in seq_along(phats)) {
    phats[[i]] <- t(apply(phats[[i]], 1, stats::quantile,
                          probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  }
  lapply(phats, as.data.frame)
}

#' @export
print.iprobitPredict <- function(x, rows = 10, digits = 3, ...) {
  if (x$type == "train") {
    cat("Training error:", paste0(decimal_place(x$error.rate, digits), "%\n"))
    cat("Brier score   :", decimal_place(x$brier.score, digits), "\n")
  } else if (!is.nan(x$error.rate)) {
    cat("Test error :", paste0(decimal_place(x$error.rate, digits), "%\n"))
    cat("Brier score:", decimal_place(x$brier.score, digits), "\n")
  } else {
    cat("Test data not provided.\n")
  }
  rows <- min(nrow(x$prob), rows)

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
  tab <- decimal_place(x$p[seq_len(rows), ], digits)
  print(tab)
  if (nrow(x$p) > rows) cat("# ... with", nrow(x$p) - rows, "more rows")
}

#' @export
print.iprobitPredict_quant <- function(x, rows = 5, digits = 3, ...) {
  rows.act <- nrow(x$prob[[1]])
  rows <- min(rows.act, rows)
  y.levels <- names(x$prob)

  brier <- decimal_place(x$brier.score, digits)
  if (x$type == "train") {
    error <- decimal_place(x$error, digits)
    tab <- rbind("Training error (%)" = error, "Brier score" = brier)
    print(as.data.frame(tab))
  } else if (all(!is.nan(x$error.rate))) {
    error <- decimal_place(x$error, digits)
    tab <- rbind("Test error (%)" = error, "Brier score" = brier)
    print(as.data.frame(tab))
  } else {
    cat("Test data not provided.\n")
  }

  for (j in seq_along(x$prob)) {
    cat("\nPredicted probabilities for Class =", y.levels[j], "\n")
    tab <- decimal_place(x$prob[[j]][seq_len(rows), ], digits)
    print(tab)
    if (rows.act > rows) cat("# ... with", rows.act - rows, "more rows\n")
  }
}
