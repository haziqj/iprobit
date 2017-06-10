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

  kernel.used <- factor(get_kernel(x))
  kernels <- levels(kernel.used)
  kernels <- gsub("FBM,", "Fractional Brownian motion with Hurst coef. ", kernels)
  x.var.list <- rep(list(NULL), length(kernels))
  x.var <- x$ipriorKernel$model$xname
  for (i in seq_along(x.var)) {
    x.var.list[[as.numeric(kernel.used[i])]] <-
      c(x.var.list[[as.numeric(kernel.used[i])]], x.var[i])
  }
  x.var.list <- lapply(x.var.list, function(x) paste0(x, collapse = ", "))
  x.var.list <- mapply(FUN = function(x, y) paste0(x, " (", y, ")"),
                       kernels, x.var.list)

  res <- list(call = x$call, kernel.used = x.var.list, tab = tab,
              maxit = x$maxit, niter = x$niter, stop.crit = x$stop.crit,
              lb = x$lower.bound, classes = x$y.levels)
  class(res) <- "iprobitSummary"
  res
}


#' @export
print.iprobitSummary <- function(x) {
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
  cat("\nVariational lower bound:", x$lb[x$niter])
  cat("\n\n")
}
