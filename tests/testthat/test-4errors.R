context("Errors")

test_that("Fit factors only", {

  dat <- gen_mixture(n = 10)
  expect_error(iprobit(as.numeric(dat$y), dat$X))

})
