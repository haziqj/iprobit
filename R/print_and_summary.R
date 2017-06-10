#' @export
print.iprobitMod <- function(x) {
  theta <- coef(x)

  cat("Lower bound value = ", x$lower.bound[x$niter], "\n")
  cat("Iterations = ", x$niter, "\n\n")
  print(round(theta, 5))
}

#' @export
summary.iprobitMod <- function(x) {
  if (is.iprobitMod_bin(x)) {
    theta <- coef(x)
    se <- x$se
  }
  if (is.iprobitMod_mult(x)) {
    tmp <- get_coef_se_mult(x)
    theta <- tmp$theta
    se <- tmp$se
  }

  tab <- cbind(
    Mean    = round(theta, digits = 4),
    S.D.    = round(se, digits = 4),
    "2.5%"  = round(theta - 1.96 * se, digits = 4),
    "97.5%" = round(theta + 1.96 * se, digits = 4)
  )

  kernel.used <- get_kernel(x)
  if (all.same(kernel.used)) kernel.used <- kernel.used[1]
  else kernel.used <- paste0(kernel.used, collapse = ", ")

  res <- list(call = x$call, kernel.used = kernel.used, tab = tab,
              maxit = x$maxit, niter = x$niter, stop.crit = x$stop.crit,
              lb = x$lower.bound)
  class(res) <- "iprobitSummary"
  res
}


#' @export
print.iprobitSummary <- function(x) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nRKHS used: ")
  cat(x$kernel.used, "\n\n")

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
