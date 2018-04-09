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

expand_lambda <- function(x, intr, intr.3plus = NULL) {
  # Helper function to expand lambda or lambda.sq (scale parameters) according
  # to any interactions specification.
  #
  # Args: lambda
  #
  # Returns: Expanded lambda.
  if (is.vector(x)) {
    # Binary models ------------------------------------------------------------
    return(iprior::.expand_Hl_and_lambda(x, x, intr, intr.3plus)$lambda)
  } else {
    # Multinomial models -------------------------------------------------------
    res <- NULL
    m <- ncol(x)
    for (j in seq_len(m)) {
      res[[j]] <- iprior::.expand_Hl_and_lambda(x[, j], x[, j], intr,
                                                intr.3plus)$lambda
    }
    return(matrix(unlist(res), ncol = m))
  }
}

# get_Hlam <- function(object, theta, theta.is.lambda = FALSE) {
#   # Obtains the kernel matrix Hlam.
#   #
#   # Args:
#   #
#   # Returns: For binary models, this calculate Hlam. For multinomial models,
#   # this calculates Hlam for every class---so a list is returned.
#   if (is.iprobit_bin(object)) {
#     if (is.matrix(theta)) theta <- theta[, 1]
#     return(iprior::.get_Hlam(object, theta, theta.is.lambda))
#   } else {
#     res <- NULL
#     m <- get_m(object)
#     for (j in seq_len(m)) {
#       res[[j]] <- iprior::.get_Hlam(object, as.numeric(theta[, j]),
#                                     theta.is.lambda)
#     }
#     return(res)
#   }
# }

get_Hlam <- function(object, theta, theta.is.lambda = FALSE) {
  # Obtains the kernel matrix Hlam.
  #
  # Args:
  #
  # Returns: Hlam
  if (is.matrix(theta)) theta <- theta[, 1]
  return(iprior::.get_Hlam(object, theta, theta.is.lambda))
}

get_Htildelam <- function(object, theta, xstar, theta.is.lambda = FALSE) {
  if (is.matrix(theta)) theta <- theta[, 1]
  return(iprior::.get_Htildelam(object, theta, xstar, theta.is.lambda))
}

get_Hlamsq <- function() {
  # Calculate Hlamsq for closed-form VB.
  #
  # Args: mod is an ipriorKernel object defined in parent environment, along
  # with lambda, lambdasq, Psql, Hsql, H2l, ind1, ind2, p, and no.int.
  #
  # Returns: Hlamsq mat for binary models, and list of Hlamsq matrices for
  # multinomial models.
  environment(Hlam_two_way_index) <- environment()
  q <- p + no.int
  if (is.matrix(lambda)) lambda <- lambda[, 1]
  if (is.matrix(lambdasq)) lambdasq <- lambdasq[, 1]

  # Calculate square terms of Hlamsq -----------------------------------------
  if (is.null(Hsql))
    square.terms <- Reduce("+", mapply("*", Psql[1:q], lambdasq[1:q],
                                       SIMPLIFY = FALSE))
  else
    square.terms <- Reduce("+", mapply("*", Hsql[1:q], lambdasq[1:q],
                                       SIMPLIFY = FALSE))

  # Calculate two-way terms of Hlamsq ----------------------------------------
  if (is.null(ind1) && is.null(ind2)) {
    two.way.terms <- 0
  } else {
    lambda.two.way <- Hlam_two_way_index(lambda, lambdasq)
    two.way.terms <- Reduce("+", mapply("*", H2l, lambda.two.way,
                                        SIMPLIFY = FALSE))
  }

  return(square.terms + two.way.terms)
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

# Helper function to determine when to stop the while loop for the VB-EM
# algorithm.
# If niter == 0 return TRUE because must complete 1 iteration.
# If niter == 1 then stop if maxit == 1, otherwise continue (nothing to compare).
# If niter > 1 then just check whether maxit reached or stop.crit reached.
loop_logical <- function() {
  lb.diff <- lb[niter] - lb[niter - 1]  # will also stop when LB becomes smaller
  crit1 <- (niter != maxit)
  crit2 <- (lb.diff > stop.crit)
  if (niter == 0) return(TRUE)
  else if (niter == 1) return(crit1)
  else return(crit1 & crit2)
}

# iprobitSE <- function(y, eta, thing1 = NULL, thing0 = NULL) {
#   if (is.null(thing1) | is.null(thing0)) {
#     thing1 <- exp(  # phi(eta) / Phi(eta)
#       dnorm(eta[y == 1], log = TRUE) - pnorm(eta[y == 1], log.p = TRUE)
#     )
#     thing0 <- -exp(  # -1 * {phi(eta) / Phi(-eta)}
#       dnorm(eta[y == 0], log = TRUE) - pnorm(-eta[y == 0], log.p = TRUE)
#     )
#   }
#
#   # Posterior variance of ystar ------------------------------------------------
#   var.ystar <- rep(NA, length(y))
#   # 1 - eta * phi(eta) / Phi(eta) - (phi(eta) / Phi(eta)) ^ 2
#   var.ystar[y == 1] <- 1 - eta[y == 1] * thing1 + (thing1 ^ 2)
#   # 1 - eta * (-1) * {phi(eta) / Phi(-eta)} - (phi(eta) / Phi(-eta)) ^ 2
#   var.ystar[y == 0] <- 1 - eta[y == 0] * thing0 + (thing0 ^ 2)
#   sqrt(var.ystar)
# }

# This is wrong
# EprodPhiZ_1 <- function(mu, sigma) {
#   Sigma <- outer(sigma, sigma, FUN = "*") + diag(1, length(mu))
#   L.inv <- solve(t(chol(Sigma)))
#   e <- L.inv %*% mu
#   prod(pnorm(e))
# }
#
# EprodPhiZ_2 <- function(mu, sigma) {
#   Sigma <- outer(sigma, sigma, FUN = "*") + diag(1, length(mu))
#   pmvnorm(upper = mu, sigma = Sigma)
# }

EprodPhiZ <- function(mu, sigma = rep(1, length(mu)), j, log = FALSE) {
  indx <- seq_along(mu)[-j]
  res <- integrate(
    function(z) {
      tmp <- 0
      for (k in indx) tmp <- tmp + pnorm(
        sigma[j] * z / sigma[k] + (mu[j] - mu[k]) / sigma[k],
        log.p = TRUE
      )
      exp(tmp) * dnorm(z)
    },
    lower = -Inf, upper = Inf
  )$value
  ifelse(isTRUE(log), log(res), res)
}

# EprodPhiZ <- function(mu, sigma) {
#   res1 <- EprodPhiZ_1(mu = mu, sigma = sigma)
#   res2 <- EprodPhiZ_2(mu = mu, sigma = sigma)
#   res3 <- EprodPhiZ_3(mu = mu, sigma = sigma)
#
#   list(res1, res2, res3)
# }
#
# k <- 20
# microbenchmark::microbenchmark(
#   method1 = EprodPhiZ_1(rnorm(k), rep(1, k)),
#   method2 = EprodPhiZ_2(rnorm(k), rep(1, k)),
#   method3 = EprodPhiZ_3(rnorm(k), rep(1, k))
# )
#
#
# EprodPhiZ(rnorm(4), rep(1, 4))
