context("Multinomial models")

test_that("Fitting multinomial models (IIA)", {

  m <- 4
	dat <- gen_mixture(n = 10, m = m)
	mod1 <- iprobit_mult(dat$y, dat$X, silent = TRUE, maxit = 5)
	expect_s3_class(mod1, "iprobitMod")
	expect_s3_class(mod1, "iprobitMod_mult")

	# Common intercept
	mod2 <- iprobit_mult(dat$y, dat$X, silent = TRUE, maxit = 5,
	                     common.intercept = TRUE)
	expect_equal(length(mod2$alpha), 1)

	# Common RKHS scale parameter
	mod3 <- iprobit_mult(dat$y, dat$X, silent = TRUE, maxit = 5,
	                     common.RKHS.scale = TRUE)
	expect_equal(length(mod3$lambda), 1)

	# Common intercept and RKHS scale parameter
	mod4 <- iprobit_mult(dat$y, dat$X, silent = TRUE, maxit = 5,
	                     common.intercept = TRUE, common.RKHS.scale = TRUE)
	expect_equal(length(mod4$alpha), 1)
	expect_equal(length(mod4$lambda), 1)

})

# test_that("Print and summary (IIA)", {
#
#   dat <- gen_mixture(n = 10)
#   mod <- iprobit_bin(dat$y, dat$X, silent = TRUE, maxit = 5)
#   expect_that(print(mod), prints_text("Lower bound value ="))
#
# })