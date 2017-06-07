context("Binary models")

test_that("Fitting binary models", {

	dat <- gen_mixture(n = 10)
	mod <- iprobit_bin(dat$y, dat$X, silent = TRUE, maxit = 5)
	expect_s3_class(mod, "iprobitMod")
	expect_s3_class(mod, "iprobitMod_bin")

})

test_that("Print and summary", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit_bin(dat$y, dat$X, silent = TRUE, maxit = 5)
  expect_that(print(mod), prints_text("Lower bound value ="))

})