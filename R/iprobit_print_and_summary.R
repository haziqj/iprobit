#' @export
print.iprobitMod <- function(x) {
  Intercept <- x$alpha
  if (length(Intercept) > 1)
    names(Intercept) <- paste0("Intercept[", seq_along(Intercept), "]")
  else
    names(Intercept) <- "Intercept"
  lambda <- x$lambda
  if (length(lambda) > 1)
    names(lambda) <- paste0("lambda[", seq_along(lambda), "]")
  else
    names(lambda) <- "lambda"
  theta <- c(Intercept, lambda)

  cat("Lower bound value = ", x$lower.bound[x$niter], "\n")
  cat("Iterations = ", x$niter, "\n\n")
  print(round(theta, 5))
}


#'
#' #' @export
#' ipriorProbitPrintAndSummary <- function(x) {
#'   y.hat <- fitted.ipriorProbit(x)$y
#'   y <- as.factor(x$y); levels(y) <- x$y.levels
#'   train.error.rate <- format(round(mean(y.hat != y) * 100, 2))
#'
#'   # Calculate 95% credibility interval for error rate --------------------------
#'   y.hat.upper <- fitted(x, "upper")$y
#'   train.error.rate.upper <- format(round(mean(y.hat.upper != x$y) * 100, 2))
#'   y.hat.lower <- fitted(x, "lower")$y
#'   train.error.rate.lower <- format(round(mean(y.hat.lower != x$y) * 100, 2))
#'
#'   list(
#'     train.error.rate = train.error.rate,
#'     error.band = c(train.error.rate.lower, train.error.rate.upper),
#'     lb = as.numeric(logLik(x))
#'   )
#' }
#'
#' #' @export
#' print.ipriorProbit <- function(x, newdata = NULL, testdata = NULL) {
#'   tmp <- ipriorProbitPrintAndSummary(x)
#'   train.error.rate <- tmp$train.error.rate
#'   cat("Training error rate:", train.error.rate, "%")
#'   if (!is.null(newdata) && !is.null(testdata)) {
#'     y.hat <- predict(x, newdata)$y
#'     test.error.rate <- format(round(mean(y.hat != testdata) * 100, 2))
#'     cat("\nTest error rate:", test.error.rate, "%")
#'   }
#' }
#'
#'
#' #' @export
#' summary.ipriorProbit <- function(x) {
#'   tmp <- ipriorProbitPrintAndSummary(x)
#'   train.error.rate <- tmp$train.error.rate
#'
#'   post.mean <- c(x$alpha, x$lambda)
#'   se <- x$se  # only for lambda alpha and lambda
#'   tab <- cbind(
#'     Mean    = round(post.mean, digits = 4),
#'     S.D.    = round(se, digits = 4),
#'     "2.5%"  = round(post.mean - 1.96 * se, digits = 4),
#'     "97.5%" = round(post.mean + 1.96 * se, digits = 4)
#'   )
#'   l <- length(x$lambda)
#'   if (l == 1) lam.name <- "lambda" else lam.name <- paste0("lambda", 1:l)
#'   rownames(tab) <- c("alpha", lam.name)
#'
#'   res <- list(tab = tab, call = x$call, kernel = x$kernel, maxit = x$maxit,
#'               stop.crit = x$stop.crit, niter = x$niter, lb = tmp$lb,
#'               error = train.error.rate, error.band = x$error.band)
#'   class(res) <- "iprobitSummary"
#'   res
#' }
#'
#' #' @export
#' print.iprobitSummary <- function(x) {
#'   cat("\nCall:\n")
#'   print(x$call)
#'   cat("\nRKHS used:", x$kernel, "\n\n")
#'   print(x$tab)
#'   if (x$niter == x$maxit) {
#'     cat("\nConvergence criterion not met. ")
#'   } else {
#'     cat("\nConverged to within", x$stop.crit, "tolerance. ")
#'   }
#'   cat("No. of iterations:", x$niter)
#'   # cat("\nModel classification error rate (%):", x$error, "within [",
#'   #     x$error.band[1], ",", x$error.band[2], "]")
#'   cat("\nModel classification error rate (%):", x$error)
#'   cat("\nVariational lower bound:", x$lb)
#'   cat("\n\n")
#' }
