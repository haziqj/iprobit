################################################################################
#
#   iprobit: Binary Probit Regression with I-priors
#   Copyright (C) 2017  Haziq Jamil
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

HlamFn <- function(env = environment()) {
  # Hl (list) and lambda (vector), both must be of same length, should be
  # defined in  environment.
  res.Hlam.mat <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
  assign("Hlam.mat", res.Hlam.mat, envir = env)
}

HlamsqFn <- function(env = environment()) {
  # Hl, Hsql, (both lists) and lambda, lambda.sq (both vectors), all of which
  # must be the same length, should be defined in  environment. Further, ind1
  # and ind2 are indices of all possible two-way multiplications obtained from
  # iprior::kernL$BlockBstuff
  if (is.null(Hsql))
    square.terms <- Reduce("+", mapply("*", Psql[1:q], lambda.sq[1:q],
                                       SIMPLIFY = FALSE))
  else
    square.terms <- Reduce("+", mapply("*", Hsql[1:q], lambda.sq[1:q],
                                       SIMPLIFY = FALSE))

  if (is.null(ind1) && is.null(ind2))
    two.way.terms <- 0
  else {
    lambda.two.way <- lambda[ind1] * lambda[ind2]
    two.way.terms <-
      Reduce("+", mapply("*", H2l, lambda.two.way, SIMPLIFY = FALSE))
    }

  res.Hlam.matsq <- square.terms + two.way.terms
  assign("Hlam.matsq", res.Hlam.matsq, envir = env)
}

is.iprobitMod_bin <- function(x) {
  any(class(x) == "iprobitMod_bin")
}

is.iprobitMod_mult <- function(x) {
  any(class(x) == "iprobitMod_mult")
}

is.iprobitData <- function(x) {
  any(class(x) == "iprobitData")
}

splitKernel <- function(kernel) {
   # Helper function to split the FBMs from the Hurst coefficients, if any
     paste(lapply(strsplit(kernel, ","), function(x) x[1]))
}

splitHurst <- function(kernel) {
  # Helper function to split the FBMs from the Hurst coefficients, if any
    suppressWarnings(
        tmp <- as.numeric(paste(lapply(strsplit(kernel, ","), function(x) x[2])))
     )
  tmp
}

get_Hurst <- function(object) {
  object$ipriorKernel$model$Hurst
}

get_one.lam <- function(object) {
  object$ipriorKernel$model$one.lam
}

get_kernel <- function(object) {
  kernel.used <- object$ipriorKernel$model$kernel
  Hurst.used <- get_Hurst(object)
  for (i in seq_along(kernel.used)) {
    if (kernel.used[i] == "FBM")
      kernel.used[i] <- paste0(kernel.used[i], ",", Hurst.used[i])
  }
  kernel.used
}

get_lambda <- function(object) {
  lambda <- object$lambda
  if (is.iprobitMod_bin(object)) {
    if (length(lambda) > 1)
      names(lambda) <- paste0("lambda[", seq_along(lambda), "]")
    else
      names(lambda) <- "lambda"
  } else if (is.iprobitMod_mult(object)) {
    if (nrow(lambda) > 1)
      rownames(lambda) <- paste0("lambda[", seq_along(lambda[, 1]), ",]")
    else
      rownames(lambda) <- "lambda"
    colnames(lambda) <- paste0("Class = ", seq_along(object$y.levels))
  } else {
    stop("Input iprobitMod objects only.")
  }

  lambda
}

get_alpha <- function(object) {
  alpha <- object$alpha
  if (is.iprobitMod_bin(object)) {
    # if (length(alpha) > 1)
    #   names(alpha) <- paste0("alpha[", seq_along(alpha), "]")
    # else
    names(alpha) <- "alpha"
  } else if (is.iprobitMod_mult(object)) {
    alpha <- matrix(alpha, nrow = 1)
    rownames(alpha) <- "alpha"
    colnames(alpha) <- paste0("Class = ", seq_along(object$y.levels))
  } else {
    stop("Input iprobitMod objects only.")
  }
  alpha
}

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

all.same <- function(v) {
  # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
  all(sapply(as.list(v[-1]), FUN = function(z) identical(z, v[1])))
}

lambdaExpand_mult <- function(x = lambda, env = iprobit.env) {
  environment(.lambdaExpand) <- env
  original.lambda <- x
  lambda.tmp <- NULL
  for (j in 1:m) {
    .lambdaExpand(original.lambda[, j], env)
    lambda.tmp[[j]] <- lambda
  }
  assign("lambda", matrix(unlist(lambda.tmp), ncol = m, nrow = l), envir = env)
}

HlamFn_mult <- function(env = environment()) {
  res.Hlam.mat <- NULL
  for (j in 1:m) {
    res.Hlam.mat[[j]] <- Reduce("+", mapply("*", Hl, lambda[, j], SIMPLIFY = FALSE))
  }
  assign("Hlam.mat", res.Hlam.mat, envir = env)
}

HlamsqFn_mult <- function(env = environment()) {
  res.Hlam.matsq <- NULL
  for (j in 1:m) {
    if (is.null(Hsql))
      square.terms <- Reduce("+", mapply("*", Psql, lambda.sq[, j],
                                         SIMPLIFY = FALSE))
    else
      square.terms <- Reduce("+", mapply("*", Hsql[1:q], lambda.sq[, j],
                                         SIMPLIFY = FALSE))

    if (is.null(ind1) && is.null(ind2))
      two.way.terms <- 0
    else {
      lambda.two.way <- lambda[ind1, j] * lambda[ind2, j]
      two.way.terms <-
        Reduce("+", mapply("*", H2l, lambda.two.way, SIMPLIFY = FALSE))
    }

    res.Hlam.matsq[[j]] <- square.terms + two.way.terms
  }

  assign("Hlam.matsq", res.Hlam.matsq, envir = env)
}

# myfun <- function() {
#   # Function to test lambdaExpand_mult() and the Hlam*_mult functions
#   set.seed(123)
#   dat <- gen_mixture(n = 4, m = 3)
#   mod <- iprior::kernL(y ~ ., dat)
#   iprobit.env <- environment()
#   list2env(mod, iprobit.env)
#   list2env(BlockBstuff, iprobit.env)
#   list2env(model, iprobit.env)
#   environment(lambdaExpand_mult) <- iprobit.env
#   environment(HlamFn_mult) <- iprobit.env
#   environment(HlamsqFn_mult) <- iprobit.env
#   m <- length(y.levels)
#   lambda <- matrix(2:3, ncol = m, nrow = l)
#   lambda.sq <- lambda ^ 2
#   lambdaExpand_mult()
#   HlamFn_mult()
#   HlamsqFn_mult()
#   list(H = Hl, lambda = lambda, Hlam.mat = Hlam.mat[[1]],
#        Hlam.matsq = Hlam.matsq[[1]])
# }

as.time <- function(x) {
  # For difftime objects
  time <- as.numeric(x)
  unit <- attr(x, "units")
  structure(list(time = time, unit = unit), class = "iprobitTime")
}

print.iprobitTime <- function(x) {
  cat(x$time, x$unit)
}

