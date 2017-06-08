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

checkLevels <- function(y) {
  y <- as.factor(y)
  y.levels <- levels(y)
  y.numeric <- as.numeric(y) - 1

  list(y = y.numeric, levels = y.levels)
}

is.iprobitMod_bin <- function(x) {
  any(class(x) == "iprobitMod_bin")
}

is.iprobitMod_mult <- function(x) {
  any(class(x) == "iprobitMod_mult")
}
