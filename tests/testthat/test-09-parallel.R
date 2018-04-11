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
#                                                par.method = "lb"),
#                     silent = TRUE)
#   mod.trer <- iprobit(y ~ ., dat, control = list(restarts = 4, maxit = 4,
#                                                  par.method = "traine"),
#                     silent = TRUE)
#   mod.tser <- iprobit(y ~ ., dat, control = list(restarts = 4, maxit = 4,
#                                                  par.method = "teste"),
#                     silent = TRUE)
#   mod.trbr <- iprobit(y ~ ., dat, control = list(restarts = 4, maxit = 4,
#                                                  par.method = "trainb"),
#                     silent = TRUE)
#   mod.tsbr <- iprobit(y ~ ., dat, control = list(restarts = 4, maxit = 4,
#                                                  par.method = "testb"),
#                     silent = TRUE)
#
#   expect_s3_class(mod.lb,   "iprobitMod")
#   expect_s3_class(mod.trer, "iprobitMod")
#   expect_s3_class(mod.tser, "iprobitMod")
#   expect_s3_class(mod.trbr, "iprobitMod")
#   expect_s3_class(mod.tsbr, "iprobitMod")
#
# })
