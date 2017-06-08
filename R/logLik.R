#' Extract the variational lower bound
#'
#' @param x an object of class \code{ipriorProbit}.
#'
#' @return The variational lower bound.
#' @export
logLik.iprobitMod <- function(x) {
  lb <- x$lower.bound[!is.na(x$lower.bound)]
  lb <- lb[length(lb)]
  class(lb) <- "iprobitLowerBound"
  lb
}

#' @export
print.iprobitLowerBound <- function(x) {
  cat("Lower bound =", x)
}
