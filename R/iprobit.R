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
                            common.intercept  = FALSE,
                            common.RKHS.scale = FALSE,
                            # nystrom = FALSE, nys.seed = NULL,
                            train.samp, control = list()) {
  # Set up controls ------------------------------------------------------------
  xname <- as.character(as.list(match.call(expand.dots = FALSE))$...)
  control_ <- list(
    maxit          = 100,
    stop.crit      = 1e-5,
    silent         = FALSE,
    alpha0         = NULL,  # if NULL, parameters are initialised in VB
    lambda0        = NULL,  # routine
    w0             = NULL,
    restarts       = 0,
    restart.method = c("lb", "error", "brier")
  )
  control <- iprior::.update_control(control, control_)
  list2env(control, environment())

  # Load the I-prior model -----------------------------------------------------
  if (iprior::is.ipriorKernel(y)) {
    mod <- y
  } else {
    mod <- iprior::kernL(y, ..., kernel = kernel, interactions = interactions,
                         # nystrom = nystrom, nys.seed = nys.seed,
                         train.samp = train.samp)
  }

  # Checks ---------------------------------------------------------------------
  if (!iprior::is.iprobit(mod)) stop("y values must be factors.")
  if (mod$no.int.3plus > 0)
    stop("Can't fit more than three-way interactions yet.")

  # Pass to the correct VB routine ---------------------------------------------
  mod$m <- m <- length(mod$y.levels)
  if (control$restarts > 1) {
    # res <- iprobit_parallel(ipriorKernel, con$restarts, con$restart.method, con)
    stop("Not implemented yet.")
  } else {
    if (m == 2) {
      res <- iprobit_bin(mod, maxit, stop.crit, silent, alpha0, lambda0, w0)
      class(res) <- c("iprobitMod", "iprobitMod_bin")
      res$coefficients <- c(get_alpha(res), get_lambda(res))
    } else {
      # res <- iprobit_mult(ipriorKernel, maxit, stop.crit, silent, alpha0,
      #                     lambda0, w0, common.intercept, common.RKHS.scale)
      # res$coefficients <- rbind(get_alpha(res), get_lambda(res))
      stop("Not implemented yet.")
    }
    res$ipriorKernel <- mod
  }

  # Change the call to "iprobit" -----------------------------------------------
  res$fullcall <- cl <- match.call()
  # ynamefromcall <- as.character(cl[2])
  # check.yname <- is.null(ipriorKernel$model$yname)
  # if (check.yname) model$yname <- ynamefromcall
  # cl[[1L]] <- as.name("iprobit")
  # m <- match(c("control"), names(cl), 0L)
  # if (any(m > 0)) cl <- cl[-m]
  res$call <- cl
  res$formula <- mod$formula

  # Include these also in the iprobitMod object --------------------------------
  res$control <- control

  res
}

#' @export
iprobit.formula <- function(formula, data = parent.frame(), kernel = "Canonical",
                            silent = FALSE, one.lam = FALSE, parsm = TRUE,
                            control = list(), ...) {
  # Pass to iprobit default ----------------------------------------------------
  ipriorKernel <- iprior::.kernL(formula, data, model = list(kernel = kernel,
                                                            one.lam = one.lam,
                                                            parsm = parsm))
  est <- iprobit.default(y = ipriorKernel, control = control, silent = silent)

  # Changing the call to simply iprobit ----------------------------------------
  cl <- match.call()
  est$fullcall <- cl
  cl[[1L]] <- as.name("iprobit")
  m <- match(c("formula", "data"), names(cl), 0L)
  cl <- cl[c(1L, m)]
  est$call <- cl
  names(est$call)[2] <- "formula"
  est$formula <- formula
  # est$terms <- class(est) <- "ipriorMod"

  est
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
