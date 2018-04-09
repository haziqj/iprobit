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



#' Extract the variational lower bound
#'
#' @param object An object of class \code{ipriorProbit}.
#' @param ... This is not used here.
#'
#' @return The variational lower bound.
#' @export
logLik.iprobitMod <- function(object, theta = NULL, silent = TRUE,
                              stop.crit = 1e-3, alpha = get_alpha(object), ...) {
  if (is.null(theta)) {
    lb <- object$lower.bound[!is.na(object$lower.bound)]
    lb <- lb[length(lb)]
    class(lb) <- "iprobitLowerBound"
    return(lb)
  } else {
    my.con <- list(
      theta0 = theta,
      alpha0 = alpha,
      stop.crit = stop.crit,
      silent = silent,
      w.only = TRUE
    )
    res <- iprobit(object$ipriorKernel, control = my.con)$lower.bound
    return(res[length(res)])
  }
}

