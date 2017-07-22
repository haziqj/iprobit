#' @export
iprobit <- function(...) {
  UseMethod("iprobit")
}

#' @export
iprobit.default <- function(y, ..., kernel = "Canonical", silent = FALSE,
                            interactions = NULL, parsm = TRUE, control = list()) {
  # Set up controls ------------------------------------------------------------
  xname <- as.character(as.list(match.call(expand.dots = FALSE))$...)
  con <- list(
    maxit             = 100,
    stop.crit         = 1e-5,
    silent            = FALSE,
    alpha0            = NULL,  # if NULL, parameters are initialised in VB
    lambda0           = NULL,  # routine
    w0                = NULL,
    common.intercept  = FALSE,
    common.RKHS.scale = FALSE,
    Nystrom           = FALSE,
    Nys.seed          = NULL
  )
  con_names <- names(con)
  con[(control_names <- names(control))] <- control
  if (length(noNms <- control_names[!control_names %in% con_names])) {
    warning("Unknown names in control options: ", paste(noNms, collapse = ", "),
            call. = FALSE)
  }
  silent_ <- silent
  list2env(con, environment())
  silent <- any(isTRUE(silent), isTRUE(silent_))

  # Pass to kernel loader and then appropriate VB routine ----------------------
  if (iprior::is.ipriorKernel(y)) {
    ipriorKernel <- y
  } else {
    ipriorKernel <- iprior::kernL(y, ...,
                                  model = list(kernel = kernel, parsm = parsm,
                                               interactions = interactions,
                                               xname = xname))
  }

  # Take samples and re-order for Nystrom --------------------------------------
  n <- ipriorKernel$n
  if (as.numeric(Nystrom) == n) Nystrom <- FALSE
  if (as.numeric(Nystrom) > 0) {
    if (!is.null(Nys.seed)) set.seed(Nys.seed)
    Nys.samp <- sample(seq_len(n), size = n, replace = FALSE)
    # in the future, should replace this simple sampling with a more
    # sophisticated procedure, preferably sampling based on each category
    ipriorKernel <- .reorder_ipriorKernel(ipriorKernel, Nys.samp)
    ipriorKernel$Nystrom <- list(m = Nystrom, Nys.samp = Nys.samp, Nys.seed = Nys.seed)
  }

  # Checks ---------------------------------------------------------------------
  if (!isTRUE(ipriorKernel$model$probit)) stop("y values must be factors.")
  if (ipriorKernel$r > 0) stop("Can't fit higher order terms yet.")
  if (ipriorKernel$no.int.3plus > 0)
    stop("Can't fit more than three-way interactions yet.")

  # Pass to the correct VB routine ---------------------------------------------
  ipriorKernel$m <- m <- length(ipriorKernel$y.levels)
  if (m == 2) {
    res <- iprobit_bin(ipriorKernel, maxit, stop.crit, silent, alpha0, lambda0, w0)
    param <- c(get_alpha(res), get_lambda(res))
  } else {
    res <- iprobit_mult(ipriorKernel, maxit, stop.crit, silent, alpha0, lambda0,
                        w0, common.intercept, common.RKHS.scale)
    param <- rbind(get_alpha(res), get_lambda(res))
  }

  # Change the call to "iprobit" -----------------------------------------------
  res$fullcall <- cl <- match.call()
  ynamefromcall <- as.character(cl[2])
  check.yname <- is.null(ipriorKernel$model$yname)
  if (check.yname) model$yname <- ynamefromcall
  cl[[1L]] <- as.name("iprobit")
  m <- match(c("control"), names(cl), 0L)
  if (any(m > 0)) cl <- cl[-m]
  res$call <- cl
  res$formula <- formula(ipriorKernel$call)

  # Include these also in the ipriorMod object ---------------------------------
  res$control      <- con
  res$coefficients <- param

  res
}

#' @export
iprobit.formula <- function(formula, data = parent.frame(), kernel = "Canonical",
                            silent = FALSE, one.lam = FALSE, parsm = TRUE,
                            control = list(), ...) {
  # Pass to iprobit default ----------------------------------------------------
  ipriorKernel <- iprior::kernL(formula, data, model = list(kernel = kernel,
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
