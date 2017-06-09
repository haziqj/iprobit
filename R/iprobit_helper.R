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

iprobitSE <- function(y, eta, thing1 = NULL, thing0 = NULL) {
  if (is.null(thing1) | is.null(thing0)) {
    thing1 <- exp(  # phi(eta) / Phi(eta)
      dnorm(eta[y == 1], log = TRUE) - pnorm(eta[y == 1], log.p = TRUE)
    )
    thing0 <- -exp(  # -1 * {phi(eta) / Phi(-eta)}
      dnorm(eta[y == 0], log = TRUE) - pnorm(-eta[y == 0], log.p = TRUE)
    )
  }

  # Posterior variance of ystar ------------------------------------------------
  var.ystar <- rep(NA, length(y))
  # 1 - eta * phi(eta) / Phi(eta) - (phi(eta) / Phi(eta)) ^ 2
  var.ystar[y == 1] <- 1 - eta[y == 1] * thing1 + (thing1 ^ 2)
  # 1 - eta * (-1) * {phi(eta) / Phi(-eta)} - (phi(eta) / Phi(-eta)) ^ 2
  var.ystar[y == 0] <- 1 - eta[y == 0] * thing0 + (thing0 ^ 2)
  sqrt(var.ystar)
}

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

EprodPhiZ <- function(mu, sigma = rep(1, length(mu)), log = FALSE) {
  res <- integrate(
    function(z) {
      tmp <- 0
      for (i in seq_len(length(mu)))
        tmp <- tmp + pnorm(z + mu[i], log.p = TRUE)
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