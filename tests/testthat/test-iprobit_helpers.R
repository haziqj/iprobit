context("iprobit helper functions")

test_that("lambda_expand_bin()", {

  lambda <- 1:3
  intr <- matrix(c(2, 3))
  res1 <- expand_lambda_bin(lambda, NULL)
  res2 <- expand_lambda_bin(lambda, intr)

})
