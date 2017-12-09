# context("Utilities and helper functions")
#
# # loop_logical() ---------------------------------------------------------------
#
# test_that("Stops when lower bound becomes smaller", {
#
#   lb <- seq(0, -100, length = 100)
#   niter <- 50
#   maxit <- 100
#   stop.crit <- 1e-5
#   environment(loop_logical) <- environment()
#   expect_false(loop_logical())
#
# })
#
# test_that("Stops when maxit reached", {
#
#   maxit <- 10
#   lb <- rep(NA, maxit)
#   niter <- 10
#   stop.crit <- 1e-5
#   environment(loop_logical) <- environment()
#   expect_false(loop_logical())
#
# })
#
# test_that("Stops after 1 iteration when maxit = 1", {
#
#   maxit <- 1
#   stop.crit <- 1e-5
#   environment(loop_logical) <- environment()
#
#   # At the beginning
#   lb <- NA
#   niter <- 0
#   expect_true(loop_logical())
#
#   # After one iteration
#   lb <- rnorm(1)
#   niter <- 1
#   expect_false(loop_logical())
#
# })
#
# test_that("Stops when stopping criterion met", {
#
#   maxit <- 1000
#   stop.crit <- 1e-5
#   environment(loop_logical) <- environment()
#   niter <- 2
#   lb <- c(pi, pi + 1e-7)
#   expect_false(loop_logical())
#
# })
