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
print.iprobitMod <- function(x, dp = 2, ...) {
  theta <- coef(x)

  cat("Training error rate:", decimal_place(x$fitted.values$train.error, dp), "%\n")
  # cat("Brier score:", decimal_place(x$brier.score, dp), "\n")
  cat("Lower bound value:", x$lower.bound[x$niter], "\n")
  # cat("Iterations = ", x$niter, "\n")
  cat("\n")
  print(round(theta, 5))
}

#' @export
summary.iprobitMod <- function(object, ...) {
  if (is.iprobitMod_bin(object)) {
    theta <- coef(object)
    se <- c(object$se.alpha, object$se.lambda)
  }
  if (is.iprobitMod_mult(object)) {
    tmp <- get_coef_se_mult(object)
    theta <- tmp$theta
    se <- tmp$se
  }

  tab <- cbind(
    Mean    = round(theta, digits = 4),
    S.D.    = round(se, digits = 4),
    "2.5%"  = round(theta - 1.96 * se, digits = 4),
    "97.5%" = round(theta + 1.96 * se, digits = 4)
  )

  kernel.used <- factor(get_kernel(object))
  kernels <- levels(kernel.used)
  kernels <- gsub("FBM,", "Fractional Brownian motion with Hurst coef. ", kernels)
  x.var.list <- rep(list(NULL), length(kernels))
  x.var <- object$ipriorKernel$model$xname
  for (i in seq_along(x.var)) {
    x.var.list[[as.numeric(kernel.used[i])]] <-
      c(x.var.list[[as.numeric(kernel.used[i])]], x.var[i])
  }
  x.var.list <- lapply(x.var.list, function(x) paste0(x, collapse = ", "))
  x.var.list <- mapply(FUN = function(x, y) paste0(x, " (", y, ")"),
                       kernels, x.var.list)

  res <- list(call = object$call, kernel.used = x.var.list, tab = tab,
              maxit = object$maxit, niter = object$niter,
              stop.crit = object$stop.crit, lb = object$lower.bound,
              classes = object$y.levels, Nystrom = object$ipriorKernel$Nystrom,
              Nystrom.check = isNystrom(object$ipriorKernel),
              train.error = object$fitted.values$train.error,
              brier.score = object$fitted.values$brier.score)
  class(res) <- "iprobitSummary"
  res
}

#' @export
print.iprobitSummary <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nClasses: ")
  cat(paste0(x$classes, collapse = ", "), "\n")

  cat("\nRKHS used:\n")
  for (i in seq_along(x$kernel.used))
    cat(x$kernel.used[[i]], "\n")

  cat("\nParameter estimates:\n")
  print(x$tab)

  if (x$niter == x$maxit) {
    cat("\nConvergence criterion not met. ")
  } else {
    cat("\nConverged to within", x$stop.crit, "tolerance. ")
  }
  cat("No. of iterations:", x$niter)

  if (isTRUE(x$Nystrom.check)) {
    cat("\nNystrom approximation used (with", x$Nystrom$m, "random subsamples)")
  }

  cat("\nVariational lower bound:", x$lb[x$niter], "\n")
  cat("Training error rate:", decimal_place(x$train.error),
      "%. Brier score:", decimal_place(x$brier.score), "\n")
  cat("\n")
}
