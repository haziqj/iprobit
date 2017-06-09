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

HlamFn <- function(env) {
  # Hl (list) and lambda (vector), both must be of same length, should be
  # defined in  environment.
  res.Hlam.mat <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
  assign("Hlam.mat", res.Hlam.mat, envir = env)
}

HlamsqFn <- function(env) {
  # Hl, Hsql, (both lists) and lambda, lambda.sq (both vectors), all of which
  # must be the same length, should be defined in  environment. Further, ind1
  # and ind2 are indices of all possible two-way multiplications obtained from
  # iprior::kernL$BlockBstuff
  if (is.null(Hsql))
    square.terms <- Reduce("+", mapply("*", Psql, lambda.sq, SIMPLIFY = FALSE))
  else
    square.terms <- Reduce("+", mapply("*", Hsql, lambda.sq, SIMPLIFY = FALSE))

  if (is.null(ind1) && is.null(ind2))
    two.way.terms <- 0
  else {
    lambda.two.way <- lambda[ind1] * lambda[ind2]
    two.way.terms <-
      Reduce("+", mapply("*", H2l, lambda.two.way, SIMPLIFY = FALSE)) +
      Reduce("+", mapply("*", lapply(H2l, t), lambda.two.way, SIMPLIFY = FALSE))
  }

  res.Hlam.matsq <<- square.terms + two.way.terms
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
