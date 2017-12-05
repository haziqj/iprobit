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
print.iprobitMod <- function(x, digits = 5, ...) {
  cat("Training error rate:", round(x$fitted.values$error, 2), "%\n")
  # cat("Brier score:", decimal_place(x$brier.score, dp), "\n")
  cat("Lower bound value:", x$lower.bound[x$niter], "\n")
  # cat("Iterations = ", x$niter, "\n")
  cat("\n")
  print(round(coef(x), digits))
}

#' @export
summary.iprobitMod <- function(object, ...) {
  tab <- object$param.summ

  kernel.used <-  iprior::.kernels_for_summary(object)

  train.error <- object$fitted.values$error
  train.brier <- object$fitted.values$brier
  test.error <- test.brier <- NULL
  if (iprior::.is.ipriorKernel_cv(object$ipriorKernel)) {
    test.error <- object$test$error
    test.brier <- object$test$brier
  }

  res <- list(tab = tab, lb = as.numeric(logLik(object)),
              train.error = train.error, train.brier = train.brier,
              test.error = test.error, test.brier = test.brier,
              call = object$call, x.kern = kernel.used,
              est.method = object$est.method, est.conv = object$est.conv,
              niter = object$niter, maxit = object$control$maxit,
              time = object$time)
  class(res) <- "iprobitMod_summary"
  res
}

#' @export
print.iprobitMod_summary <- function(x, wrap = FALSE, ...) {
  cat("Call:\n")
  cl <- x$call
  if (isTRUE(wrap)) {
    cl <- capture.output(cl)
    cl <- paste0(strwrap(cl, ...), collapse = "\n  ")
    cat(cl)
    cat("\n\n")
  } else {
    print(cl)
    cat("\n")
  }
  cat("Classes: ")
  cat(paste0(x$classes, collapse = ", "), "\n")
  cat("RKHS used:\n")
  cat(x$x.kern)
  cat("\n")
  cat("Hyperparameters:\n")
  tmp <- capture.output(print(round(x$tab, 4)))
  cat(paste(gsub("NA", "  ", tmp), collapse = "\n"))
  cat("\n")
  cat("---\n")
  cat("\n")
  cat(x$est.method)
  cat(" Iterations:", paste0(x$niter, "/", x$maxit), "\n")
  cat(x$est.conv)
  cat(" Time taken: ")
  print(x$time)
  cat("\n")
  cat("Variational lower bound:", x$lb, "\n")
  cat("Training error:", paste0(x$train.error, "%. Brier score:"), x$train.brier, "\n")
  if (!is.null(x$test.error) & !is.null(x$test.brier)) {
    cat("Test error:", paste0(x$test.error, "%. Brier score:"), x$test.brier, "\n")
  }
}
