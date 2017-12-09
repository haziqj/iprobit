context("Plots")

# Binary models ----------------------------------------------------------------

set.seed(123)
n <- 50
dat <- gen_circle(n)
mod <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "fbm",
               train.samp = sort(sample(seq_len(n), size = 40)),
               control = list(silent = TRUE, maxit = 10))

test_that("iplot_fitted()", {

  expect_silent(p <- iplot_fitted(mod))

})

test_that("iplot_dec_bound()", {

  expect_silent(p <- iplot_dec_bound(mod, grid.len = 5))

})

test_that("iplot_predict()", {

  expect_silent(p <- iplot_predict(mod, grid.len = 5))
  expect_silent(p <- iplot_predict(mod, grid.len = 5, dec.bound = FALSE))
  expect_silent(p <- iplot_predict(mod, grid.len = 5, plot.test = FALSE))
  expect_silent(p <- iplot_predict(mod, grid.len = 5, dec.bound = FALSE,
                                   plot.test = FALSE))

})

test_that("iplot_lb()", {

  expect_silent(p <- iplot_lb(mod))

})

test_that("iplot_error()", {

  expect_silent(p <- iplot_error(mod, 3))
  expect_silent(p <- iplot_error(mod, plot.test = FALSE))

})

test_that("iplot_lb_and_error()", {

  expect_silent(p <- iplot_lb_and_error(mod, niter.plot = 3, plot.test = FALSE))

})

# Binary models ----------------------------------------------------------------

set.seed(123)
n <- 100
dat <- gen_mixture(n, m = 3)
mod <- iprobit(y ~ ., dat, kernel = "fbm", train.samp = sort(sample(seq_len(n),
                                                                    size = 90)),
               control = list(silent = TRUE, maxit = 4))

test_that("iplot_fitted()", {

  expect_silent(p <- iplot_fitted(mod))

})

test_that("iplot_dec_bound()", {

  expect_silent(p <- iplot_dec_bound(mod, grid.len = 5))

})

test_that("iplot_predict()", {

  expect_silent(p <- iplot_predict(mod, grid.len = 5))
  expect_silent(p <- iplot_predict(mod, grid.len = 5, dec.bound = FALSE))
  expect_silent(p <- iplot_predict(mod, grid.len = 5, plot.test = FALSE))
  expect_silent(p <- iplot_predict(mod, grid.len = 5, dec.bound = FALSE,
                                   plot.test = FALSE))

})

test_that("iplot_lb()", {

  expect_silent(p <- iplot_lb(mod))

})

test_that("iplot_error()", {

  expect_silent(p <- iplot_error(mod, 3))
  expect_silent(p <- iplot_error(mod, plot.test = FALSE))

})

test_that("iplot_lb_and_error()", {

  expect_silent(p <- iplot_lb_and_error(mod, niter.plot = 3, plot.test = FALSE))

})

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
