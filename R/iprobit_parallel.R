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
iprobit_parallel <- function(ipriorKernel, starts,
                             method = c("lb", "error", "brier"),
                             control = list()) {
  method <- match.arg(method, c("lb", "error", "brier"))
  no.cores <- parallel::detectCores()
  if (!isTRUE(control$silent))
    cat("Performing", starts, "random restarts on", no.cores, "cores.\n")
  if (!isTRUE(control$silent)) {
    snow.options.list <- list(progress = function(i) setTxtProgressBar(pb, i))
    pb <- txtProgressBar(min = 0, max = starts, style = 1)
  } else {
    snow.options.list <- list()
  }

  # The multithreading bit -----------------------------------------------------
  cl <- parallel::makeCluster(no.cores)
  doSNOW::registerDoSNOW(cl)
  res <- foreach::`%dopar%`(
    foreach::foreach(
      i = 1:starts,
      .packages = c("iprior", "iprobit"),
      .options.snow = snow.options.list
    ), {
      new.control <- control
      new.control$maxit <- 3
      new.control$restarts <- 0
      iprobit(ipriorKernel, control = new.control, silent = TRUE)
    }
  )
  if (!isTRUE(control$silent)) close(pb)
  parallel::stopCluster(cl)

  # Find best starting value ---------------------------------------------------
  best.run <- c(
    lb    = ipar_compare_lb(res)[1],
    error = ipar_compare_error(res)[1],
    brier = ipar_compare_brier(res)[1]
  )
  if (!isTRUE(control$silent)) print(best.run)
  if (method == "lb") best.run <- best.run[1]
  if (method == "error") best.run <- best.run[2]
  if (method == "brier") best.run <- best.run[3]

  # Continue updating the best model -------------------------------------------
  mod <- res[[best.run]]
  mod$ipriorKernel <- ipriorKernel
  mod <- iprobit(mod, maxit = control$maxit - 3, silent = control$silent)
  mod
}

ipar_compare_lb <- function(x) {
  res <- sapply(x, function(z) as.numeric(logLik(z)))
  which(res == max(res))
}

ipar_compare_error <- function(x) {
  res <- sapply(x, function(z) get_error_rate(z))
  which(res == max(res))
}

ipar_compare_brier <- function(x) {
  res <- sapply(x, function(z) get_brier_score(z))
  which(res == max(res))
}
