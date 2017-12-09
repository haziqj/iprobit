# context("Multiple random starts")
#
# test_that("It works in general (silent works too)", {
#
#   dat <- gen_mixture(10)
#   mod <- iprobit(y ~ ., dat, control = list(restarts = 8, maxit = 4),
#                  silent = TRUE)
#
#   expect_s3_class(mod, "iprobitMod")
#
# })
#
# test_that("Other methods work", {
#
#   dat <- gen_mixture(10)
#   mod.lb <- iprobit(y ~ ., dat, control = list(restarts = 4, maxit = 4,
#                                                restart.method = "lb"),
#                     silent = TRUE)
#   mod.er <- iprobit(y ~ ., dat, control = list(restarts = 4, maxit = 4,
#                                                restart.method = "error"),
#                     silent = TRUE)
#   mod.br <- iprobit(y ~ ., dat, control = list(restarts = 4, maxit = 4,
#                                                restart.method = "brier"),
#                     silent = TRUE)
#
#   expect_s3_class(mod.lb, "iprobitMod")
#   expect_s3_class(mod.er, "iprobitMod")
#   expect_s3_class(mod.br, "iprobitMod")
#
# })
