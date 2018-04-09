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

#' Accessor functions for \code{ipriorMod} objects.
#'
#' @param object An \code{ipriorMod} object.
#'
#' @name Accessors
NULL

#' @describeIn Accessors Obtain all of the hyperparameters.
#' @export
get_hyperparam <- function(object) {
  check_and_get_iprobitMod(object)
  object$param.full
}

#' @describeIn Accessors Obtain the intercept.
#' @export
get_intercept <- function(object, by.class = FALSE) {
  res <- get_hyperparam(object)
  res <- res[grep("Intercept", rownames(res)), , drop = FALSE]
  if (!isTRUE(by.class)) {
    res <- c(res)
    if (is.common.intercept(object)) {
      res <- res[1]
      names(res) <- "Intercept"
    } else {
      names(res) <- get_names(object, "intercept")
    }
  }
  res
}

#' @export
get_alpha <- get_intercept

#' @export
get_lambda <- function(object, by.class = FALSE) {
  res <- get_hyperparam(object)
  res <- res[grep("lambda", rownames(res)), , drop = FALSE]
  if (!isTRUE(by.class)) {
    res <- c(res)
    if (is.common.intercept(object)) {
      res <- res[1]
      names(res) <- get_names(object, "lambda", FALSE)
    } else {
      names(res) <- get_names(object, "lambda", TRUE)
    }
  }
  res
}

#' @export
get_sd <- function(object) {
  setNames(object$param.summ$S.D., rownames(object$param.summ))
}

get_sd_alpha <- function(object) {
  res <- get_sd(object)
  res[grep("Intercept", names(res))]
}

get_sd_lambda <- function(object) {
  res <- get_sd(object)
  res[grep("lambda", names(res))]
}

#' @export
get_error_rate <- function(x) x$fitted.values$error

#' @export
get_error_rates <- function(x) {
  res <- x$error
  names(res) <- seq_along(res)
  res
}

#' @export
get_brier_score <- function(x) x$fitted.values$brier

#' @export
get_brier_scores <- function(x) {
  res <- x$brier
  names(res) <- seq_along(res)
  res
}

get_m <- function(object) {
  if (is.iprobitMod(object)) object <- object$ipriorKernel
  length(object$y.levels)
}

#' @export
get_lbs <- function(x) x$lower.bound

#' @export
print.iprobitLowerBound <- function(x, ...) {
  cat("Lower bound =", x)
}

get_theta <- function(object) object$theta

get_w <- function(object) object$w

get_y <- function(object) {
  res <- factor(object$ipriorKernel$y)
  levels(res) <- object$ipriorKernel$y.levels
  res
}
