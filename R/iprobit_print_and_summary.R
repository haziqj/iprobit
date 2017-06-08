interceptAndLambdaNames <- function(x) {
  Intercept <- x$alpha
  if (length(Intercept) > 1)
    names(Intercept) <- paste0("Intercept[", seq_along(Intercept), "] ")
  else
    names(Intercept) <- "Intercept"
  lambda <- x$lambda
  if (length(lambda) > 1)
    names(lambda) <- paste0("lambda[", seq_along(lambda), "] ")
  else
    names(lambda) <- "lambda"
  c(Intercept, lambda)
}

#' @export
print.iprobitMod <- function(x) {
  theta <- interceptAndLambdaNames(x)

  cat("Lower bound value = ", x$lower.bound[x$niter], "\n")
  cat("Iterations = ", x$niter, "\n\n")
  print(round(theta, 5))
}

#' @export
summary.iprobitMod <- function(x) {
  theta <- interceptAndLambdaNames(x)
  se <- x$se

  tab <- cbind(
    Mean    = round(theta, digits = 4),
    S.D.    = round(se, digits = 4),
    "2.5%"  = round(theta - 1.96 * se, digits = 4),
    "97.5%" = round(theta + 1.96 * se, digits = 4)
  )

  res <- list(call = x$call, kernel = x$kernel, Hurst = x$Hurst, tab = tab,
              maxit = x$maxit, niter = x$niter, stop.crit = x$stop.crit,
              lb = x$lower.bound)
  class(res) <- "iprobitSummary"
  res
}

#' @export
print.iprobitSummary <- function(x) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nRKHS used:", x$kernel)
  if (x$kernel == "FBM")
    cat(" with Hurst coefficient", x$Hurst, "\n\n")
  else
    cat("\n\n")

  cat("Parameter estimates:\n")
  print(x$tab)

  if (x$niter == x$maxit) {
    cat("\nConvergence criterion not met. ")
  } else {
    cat("\nConverged to within", x$stop.crit, "tolerance. ")
  }
  cat("No. of iterations:", x$niter)
  cat("\nVariational lower bound:", x$lb[x$niter])
  cat("\n\n")
}
