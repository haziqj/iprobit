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
