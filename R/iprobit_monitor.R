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
# # NOT USED YET
# iprobit_monitor <- function(formula, data, test.data = NULL, kernel = "Canonical",
#                             one.lam = FALSE, parsm = TRUE, monitoring = 10) {
#   for (i in seq_len(monitoring)) {
#     mod <- iprobit(formula, data, kernel = kernel, one.lam = one.lam, parsm = parsm,
#                    control = list(maxit = 1))
#     p1 <- iplot_predict(mod)
#     print(p1)
#   }
# }
