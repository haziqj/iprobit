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
is.iprobitMod_bin <- function(x) inherits(x, "iprobitMod_bin")

#' @export
is.iprobitMod_mult <- function(x) inherits(x, "iprobitMod_mult")

#' @export
is.iprobitData <- function(x) inherits(x, "iprobitData")

isNystrom <- function(x) {
  if (inherits(x, "ipriorKernel_old")) {
    if (!is.list(x$Nystrom)) res <- x$Nystrom
    else res <- TRUE
  } else {
    if (!is.list(x$ipriorKernel$Nystrom)) res <- x$ipriorKernel$Nystrom
    else res <- TRUE
  }
  res
}

#' Extract the variational lower bound
#'
#' @param object An object of class \code{ipriorProbit}.
#' @param ... This is not used here.
#'
#' @return The variational lower bound.
#' @export
logLik.iprobitMod <- function(object, ...) {
  lb <- object$lower.bound[!is.na(object$lower.bound)]
  lb <- lb[length(lb)]
  class(lb) <- "iprobitLowerBound"
  lb
}

#' @export
get_lbs <- function(x) x$lower.bound

#' @export
print.iprobitLowerBound <- function(x, ...) {
  cat("Lower bound =", x)
}

get_theta <- function(object) object$theta

get_w <- function(object) object$w

get_y <- function(object) {
  res <- factor(object$ipriorKernel$y)
  levels(res) <- object$ipriorKernel$y.levels
  res
}


#' @export
get_one.lam <- function(object) {
  object$ipriorKernel$model$one.lam
}

#' @export
get_kernel <- function(object, collapse = TRUE) {
  kernel.used <- object$ipriorKernel$model$kernel
  Hurst.used <- get_Hurst(object)
  for (i in seq_along(kernel.used)) {
    if (kernel.used[i] == "FBM")
      kernel.used[i] <- paste0(kernel.used[i], ",", Hurst.used[i])
  }
  kernel.used
}

#' get_lambda <- function(object) {
#'   lambda <- object$lambda
#'   if (is.iprobitMod_bin(object)) {
#'     if (length(lambda) > 1)
#'       names(lambda) <- paste0("lambda[", seq_along(lambda), "]")
#'     else
#'       names(lambda) <- "lambda"
#'   } else if (is.iprobitMod_mult(object)) {
#'     if (nrow(lambda) > 1)
#'       rownames(lambda) <- paste0("lambda[", seq_along(lambda[, 1]), ",]")
#'     else
#'       rownames(lambda) <- "lambda"
#'     colnames(lambda) <- paste0("Class = ", seq_along(object$y.levels))
#'   } else {
#'     stop("Input iprobitMod objects only.")
#'   }
#'
#'   lambda
#' }


# get_alpha <- function(object) {
#   alpha <- object$alpha
#   if (is.iprobitMod_bin(object)) {
#     # if (length(alpha) > 1)
#     #   names(alpha) <- paste0("alpha[", seq_along(alpha), "]")
#     # else
#     names(alpha) <- "alpha"
#   } else if (is.iprobitMod_mult(object)) {
#     alpha <- matrix(alpha, nrow = 1)
#     rownames(alpha) <- "alpha"
#     colnames(alpha) <- paste0("Class = ", seq_along(object$y.levels))
#   } else {
#     stop("Input iprobitMod objects only.")
#   }
#   alpha
# }

get_alpha <- function(object) {
  param.full <- param.summ_to_param.full(object$param.summ)
  param.full[grep("alpha", names(param.full))]
}

get_lambda <- function(object) {
  param.full <- param.summ_to_param.full(object$param.summ)
  param.full[grep("lambda", names(param.full))]
}

get_sd <- function(object) {
  setNames(object$param.summ$S.D., rownames(object$param.summ))
}

get_sd_alpha <- function(object) {
  res <- get_sd(object)
  res[grep("alpha", names(res))]
}

get_sd_lambda <- function(object) {
  res <- get_sd(object)
  res[grep("lambda", names(res))]
}

#' @export
get_coef_se_mult <- function(object) {
  theta <- coef(object)
  m <- ncol(theta)
  l <- nrow(theta) - 1

  if (isTRUE(object$control$common.intercept)) {
    alpha <- theta[1, 1]
    names(alpha) <- "alpha"
    alpha.se <- object$se.alpha
  } else {
    alpha <- theta[1, ]
    names(alpha) <- paste0("alpha[", 1:m, "]")
    alpha.se <- rep(object$se.alpha, m)
  }

  if (isTRUE(object$control$common.RKHS.scale)) {
    lambda <- theta[-1, 1]
    if (length(lambda) > 1)
      names(lambda) <- paste0("lambda[", seq_along(lambda), ",]")
    else
      names(lambda) <- "lambda"
    lambda.se <- object$se.lambda[, 1]
  } else {
    lambda <- c(t(theta[-1, ]))
    lambda.names <- NULL
    for (k in 1:l) {
      lambda.names <- c(lambda.names, paste0("lambda[", k, ",", 1:m, "]"))
    }
    names(lambda) <- lambda.names
    lambda.se <- c(t(object$se.lambda))
  }

  list(theta = c(alpha, lambda), se = c(alpha.se, lambda.se))
}

#' @export
get_error_rate <- function(x) x$fitted.values$train.error

#' @export
get_error_rates <- function(x) {
  res <- x$error
  names(res) <- seq_along(res)
  res
}

#' @export
get_brier_score <- function(x) x$fitted.values$brier.score

#' @export
get_brier_scores <- function(x) {
  res <- x$brier
  names(res) <- seq_along(res)
  res
}

all.same <- function(v) {
  # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
  all(sapply(as.list(v[-1]), FUN = function(z) identical(z, v[1])))
}

# lambda expansion and Hlam.mat calculation for binary models ------------------

lambdaExpand_bin <- function(x = lambda, y = lambda.sq, env = iprobit.env) {
  environment(.lambdaExpand) <- environment()
  original.lambda <- x
  .lambdaExpand(x = original.lambda, env = environment())
  lambda.tmp <- lambda
  assign("lambda", lambda.tmp, envir = env)

  if (!is.null(y)) {
    original.lambda.sq <- y
    lambda.sq.tmp <- NULL
    .lambdaExpand(x = original.lambda.sq, env = environment())
    lambda.sq.tmp <- lambda
    assign("lambda.sq", lambda.sq.tmp, envir = env)
  }
}

expand_lambda <- function(x, intr, intr.3plus = NULL) {
  # Helper function to expand lambda or lambda.sq (scale
  # parameters) according to any interactions specification.
  #
  # Args: lambda
  #
  # Returns: Expanded lambda


  # if lambda is vector then it is binary
  iprior::.expand_Hl_and_lambda(x, x, intr, intr.3plus)$lambda
  # else if lambda is matrix then it is multinomial
}


HlamFn <- function(env = environment()) {
  # Hl (list) and lambda (vector), both must be of same length, should be
  # defined in  environment.
  res.Hlam.mat <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
  assign("Hlam.mat", res.Hlam.mat, envir = env)
}

get_Hlam <- function(object, theta, theta.is.lambda = FALSE) {
  if (is.iprobit_bin(object)) {
    return(iprior::.get_Hlam(object = object, theta = theta,
                             theta.is.lambda = theta.is.lambda))
  } else {
    stop("Not implemented yet.")
  }
}

get_Hlamsq <- function() {
  # Args: mod is an ipriorKernel object. Needs lambda, lambdasq, Psql, Hsql,
  # H2l, ind1, ind2, p, and no.int defined in the parent environment.
  environment(Hlam_two_way_index) <- environment()
  q <- p + no.int

  if (is.iprobit_bin(mod)) {
    # Calculate square terms of Hlamsq -----------------------------------------
    if (is.null(Hsql))
      square.terms <- Reduce("+", mapply("*", Psql[1:q], lambda.sq[1:q],
                                         SIMPLIFY = FALSE))
    else
      square.terms <- Reduce("+", mapply("*", Hsql[1:q], lambda.sq[1:q],
                                         SIMPLIFY = FALSE))

    # Calculate two-way terms of Hlamsq ----------------------------------------
    if (is.null(ind1) && is.null(ind2)) {
      two.way.terms <- 0
    } else {
      lambda.two.way <- Hlam_two_way_index(lambda, lambda.sq)
      two.way.terms <- Reduce("+", mapply("*", H2l, lambda.two.way,
                                          SIMPLIFY = FALSE))
    }

    return(square.terms + two.way.terms)
  } else {
    stop("Not implemented yet.")
  }
}


# HlamsqFn <- function(env = environment()) {
#   # Hl, Hsql, (both lists) and lambda, lambda.sq (both vectors), all of which
#   # must be the same length, should be defined in  environment. Further, ind1
#   # and ind2 are indices of all possible two-way multiplications obtained from
#   # iprior::.kernL$BlockBstuff
#   environment(Hlam_two_way_index) <- env
#   if (is.null(Hsql))
#     square.terms <- Reduce("+", mapply("*", Psql[1:q], lambda.sq[1:q],
#                                        SIMPLIFY = FALSE))
#   else
#     square.terms <- Reduce("+", mapply("*", Hsql[1:q], lambda.sq[1:q],
#                                        SIMPLIFY = FALSE))
#
#   if (is.null(ind1) && is.null(ind2))
#     two.way.terms <- 0
#   else {
#     lambda.two.way <- Hlam_two_way_index(lambda, lambda.sq)
#     two.way.terms <-
#       Reduce("+", mapply("*", H2l, lambda.two.way, SIMPLIFY = FALSE))
#   }
#
#   res.Hlam.matsq <- square.terms + two.way.terms
#   assign("Hlam.matsq", res.Hlam.matsq, envir = env)
# }



# lambda expansion and Hlam.mat calculation for multinomial IIA models ---------

lambdaExpand_mult <- function(x = lambda, y = lambda.sq, env = iprobit.env) {
  environment(.lambdaExpand) <- environment()
  original.lambda <- x
  lambda.tmp <- NULL
  for (j in 1:m) {
    .lambdaExpand(x = original.lambda[, j], env = environment())
    lambda.tmp[[j]] <- lambda
  }
  assign("lambda", matrix(unlist(lambda.tmp), ncol = m), envir = env)

  if (!is.null(y)) {
    original.lambda.sq <- y
    lambda.sq.tmp <- NULL
    for (j in 1:m) {
      .lambdaExpand(x = original.lambda.sq[, j], env = environment())
      lambda.sq.tmp[[j]] <- lambda
    }
    assign("lambda.sq", matrix(unlist(lambda.sq.tmp), ncol = m), envir = env)
  }
}

HlamFn_mult <- function(env = environment()) {
  res.Hlam.mat <- NULL
  for (j in 1:m) {
    res.Hlam.mat[[j]] <- Reduce("+", mapply("*", Hl, lambda[, j], SIMPLIFY = FALSE))
  }
  assign("Hlam.mat", res.Hlam.mat, envir = env)
}

HlamsqFn_mult <- function(env = environment()) {
  environment(Hlam_two_way_index) <- env
  res.Hlam.matsq <- NULL
  for (j in 1:m) {
    if (is.null(Hsql))
      square.terms <- Reduce("+", mapply("*", Psql, lambda.sq[, j],
                                         SIMPLIFY = FALSE))
    else
      square.terms <- Reduce("+", mapply("*", Hsql, lambda.sq[, j],
                                         SIMPLIFY = FALSE))

    if (is.null(ind1) && is.null(ind2))
      two.way.terms <- 0
    else {
      lambda.two.way <- Hlam_two_way_index(lambda[, j], lambda.sq[, j])
      two.way.terms <-
        Reduce("+", mapply("*", H2l, lambda.two.way, SIMPLIFY = FALSE))
    }

    res.Hlam.matsq[[j]] <- square.terms + two.way.terms
  }

  assign("Hlam.matsq", res.Hlam.matsq, envir = env)
}

Hlam_two_way_index <- function(lam, lamsq) {
  # mod <- iprior::.kernL(Species ~ . ^ 2, iris)
  # iprobit.env <- environment()
  # list2env(mod, iprobit.env)
  # list2env(BlockBstuff, iprobit.env)
  # list2env(model, iprobit.env)

  comb.ind12 <- cbind(ind1, ind2)
  comb.ind12 <- split(comb.ind12, row(comb.ind12))

  replace_ind <- function(x) {
    if (any(x > p)) {
      here <- which(x > p)
      what <- x[here] - p
      res <- c(x[-here], intr[, what])
      sort(res)
    } else {
      x
    }
  }

  lam_ind <- function(x) {
    here <- which(duplicated(x))
    if (length(here > 0)) {
      what <- x[here]
      what.not <- x[x != what]
      prod(lamsq[what]) * prod(lam[what.not])
    } else {
      prod(lam[x])
    }
  }

  lapply(lapply(comb.ind12, replace_ind), lam_ind)
}

# myfun <- function() {
#   # Function to test lambdaExpand_mult() and the Hlam*_mult functions
#   set.seed(123)
#   dat <- gen_mixture(n = 4, m = 3)
#   mod <- iprior::.kernL(y ~ . ^ 2, dat)
#   iprobit.env <- environment()
#   list2env(mod, iprobit.env)
#   list2env(BlockBstuff, iprobit.env)
#   list2env(model, iprobit.env)
#   environment(lambdaExpand_mult) <- iprobit.env
#   environment(HlamFn_mult) <- iprobit.env
#   environment(HlamsqFn_mult) <- iprobit.env
#   environment(Hlam_two_way_index) <- iprobit.env
#   m <- length(y.levels)
#   lambda <- matrix(2:3, ncol = m, nrow = l)
#   lambda.sq <- lambda ^ 2
#   lambdaExpand_mult(x = lambda, y = lambda.sq, env = iprobit.env)
#   HlamFn_mult(env = iprobit.env)
#   HlamsqFn_mult(env = iprobit.env)
#   list(H = Hl, Hsq = Hsql, lambda = lambda, lambda.sq = lambda.sq,
#        Hlam.mat = Hlam.mat[[1]], Hlam.matsq = Hlam.matsq[[1]])
#
#   # CHECK
#   # 2 3 6 (lambda)
#   # 4 9 36 (lambda.sq)
#   #
#   # H1 <- H[[1]]
#   # H2 <- H[[2]]
#   # H12 <- H[[1]] * H[[2]]
#   # H1sq <- H1 %*% H1
#   # H2sq <- H2 %*% H2
#   # H12sq <- H12 %*% H12
#   #
#   # 4 * H1sq + 9 * H2sq + 36 * H12sq + 2 * 3 * (H1 %*% H2 + H2 %*% H1) + 4 * 3 * (H1 %*% H12 + H12 %*% H1) + 9 * 2 * (H2 %*% H12 + H12 %*% H2)
# }

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







