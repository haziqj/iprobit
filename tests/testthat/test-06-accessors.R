# context("Methods for iprobitMod objects")

# Update/refit -----------------------------------------------------------------
#
# test_that("Update from formula", {
#
#   dat <- gen_mixture(10)
#   mod.original <- mod <- iprobit(y ~ ., dat, control = list(maxit = 1),
#                                  silent = TRUE)
#   mod <- iprobit(mod, maxit = 5)
#
#   expect_s3_class(mod, "iprobitMod")
#   expect_equal(mod$niter, 6)
#   expect_equal(length(mod$lower.bound), 6)
#   expect_equal(length(mod$error), 6)
#   expect_equal(length(mod$brier), 6)
#   expect_equal(mod$formula, mod.original$formula)
#
# })
#
# test_that("Update from non-formula", {
#
#   dat <- gen_mixture(10)
#   mod.original <- mod <- iprobit(dat$y, dat$X, control = list(maxit = 1),
#                                  silent = TRUE)
#   mod <- iprobit(mod, maxit = 5)
#
#   expect_s3_class(mod, "iprobitMod")
#   expect_equal(mod$niter, 6)
#   expect_equal(length(mod$lower.bound), 6)
#   expect_equal(length(mod$error), 6)
#   expect_equal(length(mod$brier), 6)
#   expect_equal(mod$formula, mod.original$formula)
#
# })
#
# test_that("Update options", {
#
#   dat <- gen_mixture(10)
#   mod.original <- mod <- iprobit(dat$y, dat$X, control = list(maxit = 1),
#                                  silent = TRUE)
#
#   expect_message(mod <- iprobit(mod))
#   expect_output(mod <- iprobit(mod, stop.crit = 1e-1, maxit = 5, silent = FALSE))
#
# })
#
# test_that("The function update()", {
#
#   dat <- gen_mixture(10)
#   mod.original <- mod <- iprobit(dat$y, dat$X, control = list(maxit = 1),
#                                  silent = TRUE)
#   update(mod, maxit = 5)
#
#   expect_s3_class(mod, "iprobitMod")
#   expect_equal(mod$niter, 6)
#   expect_equal(length(mod$lower.bound), 6)
#   expect_equal(length(mod$error), 6)
#   expect_equal(length(mod$brier), 6)
#   expect_equal(mod$formula, mod.original$formula)
#
# })
