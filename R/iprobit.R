#' @export
iprobit <- function(...) {
  UseMethod("iprobit")
}

#' @export
iprobit.default <- function(y, ..., kernel = "Canonical", silent = FALSE,
                            interactions = NULL, control = list()) {
  # Set up controls ------------------------------------------------------------
  con <- list(
    maxit             = 200,
    stop.crit         = 1e-5,
    silent            = FALSE,
    alpha0            = NULL,  # if NULL, parameters are initialised in VB
    lambda0           = NULL,  # routine
    w0                = NULL,
    common.intercept  = FALSE,
    common.RKHS.scale = FALSE
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
    ipriorKernel <- iprior::kernL(y, ..., model = list(kernel = kernel))
  }
  if (!isTRUE(ipriorKernel$model$probit)) stop("y values must be factors.")

  res <- iprobit_bin(ipriorKernel, maxit, stop.crit, silent, alpha0, lambda0, w0)
  res
}

#' @export
iprobit.formula <- function(formula, data = parent.frame(), kernel = "Canonical",
                            silent = FALSE, control = list()) {
  # Pass to iprobit default ----------------------------------------------------
  ipriorKernel <- iprior::kernL(formula, data, model = list(kernel = kernel))
  est <- iprobit.default(y = ipriorKernel, control = control)

  # Changing the call to simply iprobit ----------------------------------------
  # cl <- match.call()
  # est$fullcall <- cl
  # cl[[1L]] <- as.name("iprior")
  # m <- match(c("formula", "data"), names(cl), 0L)
  # cl <- cl[c(1L, m)]
  # est$call <- cl
  # names(est$call)[2] <- "formula"
  # est$formula <- formula
  # est$terms <- class(est) <- "ipriorMod"

  est
}
