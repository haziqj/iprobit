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
iprobit_parallel <- function(mod, method = "vb",
                             control = list(silent = FALSE, restarts = TRUE,
                                            no.cores = parallel::detectCores(),
                                            par.maxit = 10)) {

  # Set up controls ------------------------------------------------------------
  if (control$restarts == 1) control$restarts <- parallel::detectCores()
  control$no.cores <- min(parallel::detectCores(), control$restarts)
  if (!is.null(control$theta0)) {
    message("Ignoring theta0 control options with random restarts.")
  }
  control$theta0 <- NULL
  if (!isTRUE(control$silent)) {
    cat("Performing", control$restarts, "random restarts on", control$no.cores,
        "cores\n")
    snow.options.list <- list(progress = function(i) setTxtProgressBar(pb, i))
    pb <- txtProgressBar(min = 0, max = control$restarts, style = 1)
  } else {
    snow.options.list <- list()
  }
  par.method <- match.arg(control$par.method, c("lower.bound", "train.error",
                                                "train.brier", "test.error",
                                                "test.brier"))

  control$par.maxit <- control$maxit  # for now, it just does multiple restarts in parallel

  # The multithreading bit -----------------------------------------------------
  start.time <- Sys.time()
  cl <- parallel::makeCluster(control$no.cores)
  doSNOW::registerDoSNOW(cl)
  res <- foreach::`%dopar%`(
    foreach::foreach(
      i = seq_len(control$restarts),
      .packages = c("iprior", "iprobit"),
      .options.snow = snow.options.list
    ), {
      new.control          <- control
      new.control$restarts <- 0
      new.control$maxit    <- control$par.maxit
      new.control$silent   <- TRUE
      tmp <- iprobit(mod, control = new.control, method = method)
      # tmp$iprioKernel <- NULL
      tmp
    }
  )
  if (!isTRUE(control$silent)) close(pb)
  parallel::stopCluster(cl)

  # Find best starting value ---------------------------------------------------
  list2env(find_best_run(res, par.method), envir = environment())
  # best.niter <- res[[best.run]]$niter
  # best.lb    <- res[[best.run]]$lb

  # Continue updating the best model -------------------------------------------
  if (!isTRUE(control$silent)) {
    cat(paste0(par.msg, " from random starts:\n"))
    print(run.res)  # obtained from find_best_run()
    cat("Continuing on Run", best.run, "\n")
  }

  # control$restarts <- 0
  # control$alpha0   <- res[[best.run]]$alpha
  # control$theta0   <- res[[best.run]]$theta
  # control$maxit    <- control$maxit - control$par.maxit
  # res <- iprobit(mod, method = method, control = control)
  # end.time <- Sys.time()
  # time.taken <- as.time(end.time - start.time)
  #
  # # Update time taken ----------------------------------------------------------
  # res$time <- time.taken
  # res$start.time <- start.time
  # res$end.time <- end.time
  # res$niter <- res$niter + best.niter
  # res$lower.bound <- c(best.lb, res$lower.bound)

  res[[best.run]]
}

# ipar_compare_lb <- function(x) {
#   res <- sapply(x, function(z) as.numeric(logLik(z)))
#   which(res == max(res))
# }
#
# ipar_compare_error <- function(x) {
#   res <- sapply(x, function(z) get_error_rate(z))
#   which(res == min(res))
# }
#
# ipar_compare_brier <- function(x) {
#   res <- sapply(x, function(z) get_brier_score(z))
#   which(res == min(res))
# }

find_best_run <- function(res, par.method) {
  # Args: res is a list coming from the foreach output. It contains alpha,
  # theta, niter, lb, traine, trainb, teste and testb.
  par.ind <- grep(par.method, names(res[[1]]))
  tmp <- sapply(res, function(x) {
    x <- x[[par.ind]]
    x[length(x)]
  })
  if (is.list(tmp)) {  # if cannot find par.method, tmp is a list not a vector
    if (par.method == "test.error") par.method <- "train.error"
    else if (par.method == "test.brier") par.method <- "train.brier"
    else par.method <- "lower.bound"
    message("Using training results as test results not found.")
    par.ind <- grep(par.method, names(res[[1]]))
    tmp <- sapply(res, function(x) {
      x <- x[[par.ind]]
      x[length(x)]
    })
  }
  names(tmp) <- paste("Run", seq_along(tmp))

  if (par.method == "lower.bound") {
    best.run <- which(tmp == max(tmp))
    par.msg <- "Variational lower-bound"
  } else if (par.method == "train.error") {
    best.run <- which(tmp == min(tmp))
    par.msg <- "Training misclassification (percent)"
  } else if (par.method == "train.brier") {
    best.run <- which(tmp == min(tmp))
    par.msg <- "Training Brier score"
  } else if (par.method == "test.error") {
    best.run <- which(tmp == min(tmp))
    par.msg <- "Test misclassification (percent)"
  } else if (par.method == "test.brier") {
    best.run <- which(tmp == min(tmp))
    par.msg <- "Test Brier score"
  }

  if (length(best.run) > 1) best.run <- best.run[sample(seq_along(best.run), 1)]

  list(run.res = tmp, best.run = best.run, par.msg = par.msg)
}
