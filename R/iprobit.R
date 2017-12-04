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
iprobit <- function(...) UseMethod("iprobit")

#' @export
iprobit.default <- function(y, ..., kernel = "linear",  interactions = NULL,
                            est.hurst = FALSE, est.lengthscale = FALSE,
                            est.offset = FALSE,
                            common.intercept  = FALSE,
                            common.RKHS.scale = FALSE,
                            # nystrom = FALSE, nys.seed = NULL,
                            train.samp, control = list()) {
  # Load the I-prior model -----------------------------------------------------
  if (iprior::is.ipriorKernel(y)) {
    mod <- y
  } else {
    mod <- iprior::kernL(y, ..., kernel = kernel, interactions = interactions,
                         est.lambda = TRUE, lambda = lambda, est.psi = TRUE,
                         psi = 1, est.hurst = est.hurst,
                         est.lengthscale = est.lengthscale,
                         est.offset = est.offset,
                         # nystrom = nystrom, nys.seed = nys.seed,
                         train.samp = train.samp)
  }

  # Set up controls ------------------------------------------------------------
  control_ <- list(
    maxit          = 100,
    stop.crit      = 1e-5,
    silent         = FALSE,
    alpha0         = NULL,  # if NULL, parameters are
    lambda0        = NULL,  # initialised in VB
    w0             = NULL,  # routine
    restarts       = 0,
    restart.method = c("lb", "error", "brier")
  )
  control <- iprior::.update_control(control, control_)
  list2env(control, environment())

  # Checks ---------------------------------------------------------------------
  if (!iprior::.is.categorical(mod)) stop("y values must be factors.")
  if (mod$no.int.3plus > 0)
    stop("Can't fit more than three-way interactions yet.")
  mod$m <- m <- length(mod$y.levels)  # no. of classes
  est.method <- iprior::.iprior_method_checker(mod, "em")
  # >>> CHECK IF PSI = 1 <<<

  # Pass to the correct VB routine ---------------------------------------------
  if (control$restarts > 1) {
    # res <- iprobit_parallel(ipriorKernel, con$restarts, con$restart.method, con)
    stop("Not implemented yet.")
  } else {
    if (m == 2) {
      if (est.method["em.closed"]) {  # VB CLOSED-FORM
        res <- iprobit_bin(mod, maxit, stop.crit, silent, alpha0, lambda0, w0)
        res$est.method <- "Closed-form VB-EM algorithm."
      } else {
        stop("Not implemented yet.")
      }
      class(res) <- c("iprobitMod", "iprobitMod_bin")
      res$coefficients <- iprior::.reduce_theta(
        param.summ_to_param.full(res$param.summ), mod$estl
      )$theta.reduced
    } else {
      # res <- iprobit_mult(ipriorKernel, maxit, stop.crit, silent, alpha0,
      #                     lambda0, w0, common.intercept, common.RKHS.scale)
      # res$coefficients <- rbind(get_alpha(res), get_lambda(res))
      stop("Not implemented yet.")
    }
    if (res$conv == 0)
      res$est.conv <- paste0("Converged to within ", control$stop.crit,
                             " tolerance.")
    else if (res$conv == 1)
      res$est.conv <- "Convergence criterion not met."
    else
      res$est.conv <- res$message
    mod$estl$est.psi <- FALSE # <<<<<<<<<<<<<<<<<<<<<<<<<<<
    res$ipriorKernel <- mod
  }

  # Change the call to "iprobit" -----------------------------------------------
  res$call <- iprior::.fix_call_default(match.call(), "iprobit")
  res$ipriorKernel$call <- iprior::.fix_call_default(match.call(), "kernL")

  # Include these also in the iprobitMod object --------------------------------
  res$control <- control

  res
}

#' @export
iprobit.formula <- function(formula, data, kernel = "linear", one.lam = FALSE,
                            est.hurst = FALSE, est.lengthscale = FALSE,
                            est.offset = FALSE, common.intercept  = FALSE,
                            common.RKHS.scale = FALSE, lambda = 1,
                            # nystrom = FALSE, nys.seed = NULL,
                            train.samp, control = list(), ...) {
  # Simply load the kernel and pass to iprobit.default() ------------------------
  mod <- iprior::kernL(formula, data, kernel = kernel, one.lam = one.lam,
                       est.lambda = TRUE, est.hurst = est.hurst,
                       est.lengthscale = est.lengthscale,
                       est.offset = est.offset, est.psi = FALSE,
                       lambda = lambda, psi = 1,
                       # nystrom = nystrom, nys.seed = nys.seed,
                       train.samp = train.samp, ...)
  res <- iprobit.default(y = mod, control = control)
  res$call <- iprior::.fix_call_formula(match.call(), "iprobit")
  res$ipriorKernel$call <- iprior::.fix_call_formula(match.call(), "kernL")
  res
}

#' @export
iprobit.iprobitMod <- function(object, maxit = NULL, stop.crit = NULL,
                               silent = NULL, ...) {
  ipriorKernel <- object$ipriorKernel
  con          <- object$control
  con$w0       <- object$w
  con$lambda0  <- object$lambda
  con$alpha0   <- object$alpha
  con$restarts <- 0
  if (!is.null(maxit)) con$maxit <- maxit
  else {
    con$maxit <- 100
    message("Updating iprobit model with 100 iterations.")
  }
  if (!is.null(stop.crit)) con$stop.crit <- stop.crit
  if (!is.null(silent)) con$silent <- silent

  # Pass to iprobit.default ----------------------------------------------------
  res <- iprobit.default(y = ipriorKernel, control = con)

  # Update time, call, maxit, niter, lb, error, brier --------------------------
  new.time.diff <- res$end.time - res$start.time
  old.time.diff <- object$end.time - object$start.time
  res$time <- iprior::as.time(new.time.diff + old.time.diff)
  res$end.time <- object$end.time + new.time.diff
  res$call <- object$call
  res$control$maxit <- res$maxit <- res$maxit + object$maxit
  res$niter <- res$niter + object$niter
  res$lower.bound <- c(object$lower.bound, res$lower.bound)
  res$error <- c(object$error, res$error)
  res$brier <- c(object$brier, res$brier)

  res
}

#' @export
update.iprobitMod <- function(object, maxit = NULL, stop.crit = NULL,
                              silent = NULL, ...) {
  res <- iprobit.iprobitMod(object, maxit, stop.crit, silent, ...)
  assign(deparse(substitute(object)), res, envir = parent.frame())
}
