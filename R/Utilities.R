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

remove_psi <- function(object) {
  # In iprobit fitting, the model error precision psi is fixed at 1. This sets
  # psi = 1 if it was not. Most relevant for fitting through
  # iprobit.ipriorKernel().
  #
  # Args: An ipriorKernel object.
  #
  # Returns: An updated ipriorKernel object. A warning message if psi was set to
  # be estimated.
  if (isTRUE(object$estl$est.psi)) {
    warning("iprobit does not estimate model error precision. Setting psi = 1.",
            call. = FALSE, immediate. = TRUE)
    object$estl$est.psi <- FALSE
    object$thetal$theta <- object$thetal$theta[!grepl("psi", names(object$thetal$theta))]
    object$thetal$theta.omitted <- c(object$thetal$theta.omitted, "psi" = 0)
    object$thetal$theta.drop["psi"] <- TRUE
    object$thetal$n.theta <- object$thetal$n.theta - 1
  }
  object
}

#' @export
logit <- function(x, exp.logit = FALSE) {
  res <- log(x) - log(1 - x)
  if (isTRUE(exp.logit)) exp(res)
  else res
}

#' @export
expit <- function(x, log.expit = FALSE) {
  res <- -log(1 + exp(-x))
  if (isTRUE(log.expit)) res
  else exp(res)
}

is.iprobit <- function(x) {
  if (iprior::is.ipriorKernel(x)) return(x$iprobit)
  else stop("Not an ipriorKernel object.")
}

is.iprobit_bin <- function(x) {
  if (iprior::.is.categorical(x)) {
    length(x$y.levels) == 2
  } else {
    FALSE
  }
}

#' @export
is.iprobitMod <- function(x) inherits(x, "iprobitMod")

#' @export
is.iprobitMod_bin <- function(x) inherits(x, "iprobitMod_bin")

#' @export
is.iprobitMod_mult <- function(x) inherits(x, "iprobitMod_mult")

#' @export
is.iprobitData <- function(x) inherits(x, "iprobitData")

# isNystrom <- function(x) {
#   if (inherits(x, "ipriorKernel_old")) {
#     if (!is.list(x$Nystrom)) res <- x$Nystrom
#     else res <- TRUE
#   } else {
#     if (!is.list(x$ipriorKernel$Nystrom)) res <- x$ipriorKernel$Nystrom
#     else res <- TRUE
#   }
#   res
# }

is.common.intercept <- function(x) {
  check_and_get_iprobitMod(x)
  return(x$common$intercept)
}

is.common.RKHS.param <- function(x) {
  check_and_get_iprobitMod(x)
  return(x$common$RKHS.param)
}

check_and_get_iprobitMod <- function(object, assign.to.env = FALSE) {
  # Helper function to check whether object is of ipriorMod class.
  #
  # Args: An ipriorMod or ipriorKernel object; logical assign.to.env.
  #
  # Returns: Nothing - just checks. Unless assign.to.env is TRUE.
  if (is.iprobitMod(object)) {
    if (isTRUE(assign.to.env)) list2env(object$ipriorKernel, parent.frame())
  } else {
    stop("Input an iprobitMod object.", call. = FALSE)
  }
}

all.same <- function(v) {
  # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
  all(sapply(as.list(v[-1]), FUN = function(z) identical(z, v[1])))
}

hyperparam_to_theta <- function(lambda, hurst, lengthscale, offset) {
  # In order to be compatible with the iprior helper functions, need theta (the
  # unconstrained versions of the hyperparameters, excluding alpha).
  res <- NULL
  if (!missing(lambda)) {
    # if (length(lambda == 1)) res <- c(res, log(lambda))
    res <- c(res, lambda)
  }
  if (!missing(hurst)) {
    res <- c(res, qnorm(hurst))
  }
  if (!missing(lengthscale)) {
    res <- c(res, log(lengthscale))
  }
  if (!missing(offset)) {
    res <- c(res, log(offset))
  }
  res
}

theta_to_hyperparam <- function(theta, alpha, object) {
  # Convert theta to hyperparam = c(alpha, lambda, etc.)
  res <- iprior::.theta_to_collapsed_param(theta, object)
  res <- iprior::.reduce_theta(res, object$estl)$theta.reduced
  if (length(alpha) == 1) names(alpha) <- "alpha"
  # res <- res[-grep("psi", names(res))]
  c(alpha, res)
}

param.summ_to_param.full <- function(param.summ) {
  setNames(param.summ$Mean, rownames(param.summ))
}

expand_param.summ <- function(object, param.summ, theta) {
  # Args: An ipriorKernel object and param.summ table.
  #
  # Returns: An expanded param.summ table.
  alpha.ind <- grep("alpha", rownames(param.summ))
  res1 <- param.summ[alpha.ind, , drop = FALSE]
  res2 <- as.list(param.summ[-alpha.ind, , drop = FALSE])
  tmp1 <- iprior::.theta_to_collapsed_param(theta, object)
  tmp2 <- lapply(res2[-1], iprior::.expand_theta, theta.omitted = NA,
                 theta.drop = object$thetal$theta.drop)
  res2 <- as.data.frame(c(list(Mean = tmp1), tmp2))
  colnames(res2) <- colnames(res1)
  rbind(res1, res2)
}

theta_to_coef <- function(theta, object) {
  # Args: An ipriorKernel object.
  #
  # Returns: Hyperparameters of the kernel, excluding psi usually.
  param.full <- iprior::.theta_to_collapsed_param(theta, object)
  iprior::.reduce_theta(param.full, object$estl)$theta.reduced
}

theta_to_param.full <- function(theta, alpha, object) {
  # Convert theta into param full.
  #
  # Args: theta is in a table form, each column representing each class. alpha
  # is a vector of length m (no. of classes) and object is an ipriorKernel
  # object.
  #
  # Returns: "Collapsed param" but in tabular form, with each column
  # representing each class.
  res <- apply(theta, 2, iprior::.theta_to_collapsed_param, object = object)
  rownames(res) <- gsub("]", ",]", rownames(res))
  res <- rbind("Intercept" = alpha, res)
  colnames(res) <- paste0("Class = ", seq_len(ncol(res)))
  res
}

param.full_to_coef <- function(param.full, object) {
  # Convert param.full (tabular form) to coefficients (tabular form).
  #
  # Args: param.full (tabular form) and ipriorKernel object.
  #
  # Returns: Coefficients in tabular form.
  res <- apply(param.full[-1, ], 2, function(x) {
    iprior::.reduce_theta(x, est.list = object$estl)$theta.reduced
  })
  res <- rbind(param.full[1, ], res)

  tmp.names <- rownames(param.full)[param.full[, 1] %in% res[, 1]]
  rownames(res) <- tmp.names
  res
}

get_names <- function(object, names = c("intercept", "lambda", "hurst",
                                        "lengthscale", "offset"),
                      expand = TRUE) {
  m <- get_m(object)
  full.names <- names(object$thetal$theta.drop)
  lambda.count <- sum(grepl("lambda", full.names))
  hurst.count <- sum(grepl("hurst", full.names))
  lengthscale.count <- sum(grepl("lengthscale", full.names))
  offset.count <- sum(grepl("offset", full.names))

  if (isTRUE(expand)) {
    alpha.names <- paste0("Intercept[", 1:m, "]")
  } else {
    alpha.names <- "Intercept"
  }

  if (isTRUE(expand)) {
    if (lambda.count > 1) {
      lambda.ind <- expand.grid(seq_len(m), seq_len(lambda.count))
      lambda.names <- paste0("lambda[", lambda.ind[, 2], ",", lambda.ind[, 1], "]")
    } else {
      lambda.names <- paste0("lambda[", 1:m, "]")
    }
  } else {
    if (lambda.count > 1) {
      lambda.names <- paste0("lambda[", seq_len(lambda.count), ",]")
    } else {
      lambda.names <- paste0("lambda")
    }
  }

  if (hurst.count > 0) {
    if (isTRUE(expand)) {
      if (hurst.count > 1) {
        hurst.ind <- expand.grid(seq_len(m), seq_len(hurst.count))
        hurst.names <- paste0("hurst[", hurst.ind[, 2], ",", hurst.ind[, 1], "]")
      } else {
        hurst.names <- paste0("hurst[", 1:m, "]")
      }
    } else {
      if (hurst.count > 1) {
        hurst.names <- paste0("hurst[", seq_len(hurst.count), ",]")
      } else {
        hurst.names <- paste0("hurst")
      }
    }
  } else {
    hurst.names <- NULL
  }

  if (lengthscale.count > 0) {
    if (isTRUE(expand)) {
      if (lengthscale.count > 1) {
        lengthscale.ind <- expand.grid(seq_len(m), seq_len(lengthscale.count))
        lengthscale.names <- paste0("lengthscale[", lengthscale.ind[, 2], ",",
                                    lengthscale.ind[, 1], "]")
      } else {
        lengthscale.names <- paste0("lengthscale[", 1:m, "]")
      }
    } else {
      if (lengthscale.count > 1) {
        lengthscale.names <- paste0("lengthscale[", seq_len(lengthscale.count), ",]")
      } else {
        lengthscale.names <- paste0("lengthscale")
      }
    }
  } else {
    lengthscale.names <- NULL
  }

  if (offset.count > 0) {
    if (isTRUE(expand)) {
      if (offset.count > 1) {
        offset.ind <- expand.grid(seq_len(m), seq_len(offset.count))
        offset.names <- paste0("offset[", offset.ind[, 2], ",", offset.ind[, 1], "]")
      } else {
        offset.names <- paste0("offset[", 1:m, "]")
      }
    } else {
      if (offset.count > 1) {
        offset.names <- paste0("offset[", seq_len(offset.count), ",]")
      } else {
        offset.names <- paste0("offset")
      }
    }
  } else {
    offset.names <- NULL
  }

  if (isTRUE(expand)) {
    psi.names <- paste0("psi[", 1:m, "]")
  } else {
    psi.names <- "psi"
  }

  res.names <- c("intercept", "lambda", "hurst", "lengthscale", "offset", "psi")
  res <- list(alpha.names, lambda.names, hurst.names, lengthscale.names,
              offset.names, psi.names)
  res.ind <- res.names %in% names
  unlist(res[res.ind])
}


