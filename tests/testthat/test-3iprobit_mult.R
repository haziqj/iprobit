context("Multinomial models")

test_that("Fitting multinomial models (IIA)", {

  m <- 4
	dat <- gen_mixture(n = 10, m = m)
	mod1 <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
	modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE,
	                control = list(maxit = 5))
	expect_s3_class(mod1, "iprobitMod")
	expect_s3_class(mod1, "iprobitMod_mult")
	expect_s3_class(modf, "iprobitMod")
	expect_s3_class(modf, "iprobitMod_mult")

	# Common intercept
	mod2 <- iprobit(dat$y, dat$X, silent = TRUE,
	                control = list(maxit = 5, common.intercept = TRUE))
	expect_true(all.same(get_alpha(mod2)))

	# Common RKHS scale parameter
	mod3 <- iprobit(dat$y, dat$X, silent = TRUE,
	                control = list(maxit = 5, common.RKHS.scale = TRUE))
	expect_true(all.same(get_lambda(mod3)))

	# Common intercept and RKHS scale parameter
	mod4 <- iprobit(dat$y, dat$X, silent = TRUE,
	                control = list(maxit = 5,
	                               common.intercept = TRUE,
	                               common.RKHS.scale = TRUE))
	expect_true(all.same(get_alpha(mod4)))
	expect_true(all.same(get_lambda(mod4)))

})

test_that("Print (IIA)", {

  m <- 4
  dat <- gen_mixture(n = 10, m = m)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE, control = list(maxit = 5))
  expect_that(print(mod), prints_text("Lower bound value ="))
  expect_that(print(modf), prints_text("Lower bound value ="))

})

test_that("Summary", {

  m <- 4
  dat <- gen_mixture(n = 10, m = m)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE, control = list(maxit = 5))
  mod.summary <- summary(mod)
  modf.summary <- summary(modf)
  expect_s3_class(mod.summary, "iprobitSummary")
  expect_s3_class(modf.summary, "iprobitSummary")

})

test_that("Fitted", {

  m <- 4
  dat <- gen_mixture(n = 10, m = m)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE, control = list(maxit = 5))
  expect_that(fitted(mod), is_a("list"))
  expect_that(fitted(modf), is_a("list"))

})

test_that("Predict (without test error rate)", {

  m <- 4
  dat <- gen_mixture(n = 10, m = m)
  dat.test <- gen_mixture(n = 5, m = m)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE,
                  control = list(maxit = 5))

  mod.predict <- predict(mod, newdata = list(dat.test$X))
  modf.predict <- predict(modf, newdata = as.data.frame(dat.test)[, -3])

  expect_s3_class(mod.predict, "iprobitPredict")
  expect_s3_class(modf.predict, "iprobitPredict")
  expect_that(print(mod.predict), prints_text("Test data not provided."))

})

test_that("Predict (with test error rate)", {

  m <- 4
  dat <- gen_mixture(n = 10, m = m)
  dat.test <- gen_mixture(n = 5, m = m)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE,
                  control = list(maxit = 5))

  mod.predict <- predict(mod, newdata = list(dat.test$X), y = dat.test$y)
  modf.predict <- predict(modf, newdata = dat.test)

  expect_that(print(mod.predict), prints_text("Test error rate"))
  expect_that(print(modf.predict), prints_text("Test error rate"))

})

# test_that("Convergence", {
#
#   set.seed(123)
#   m <- 3
#   dat <- gen_mixture(n = 10, m = m)
#   mod <- iprobit(dat$y, dat$X, control = list(maxit = 500, silent = TRUE))
#   modf <- iprobit(y ~ ., dat, one.lam = TRUE,
#                   control = list(maxit = 500, silent = TRUE))
#   expect_equal(mod$lambda, 0.11646, tolerance = 1e-3)
#   expect_equal(modf$lambda, 0.11646, tolerance = 1e-3)
#
#   # Single lambda
#   # mod <- iprobit_bin(dat$y, dat$X, silent = TRUE, maxit = 200)
#   # > mod
#   # Lower bound value =  -4.89696
#   # Iterations =  106
#   #
#   # alpha   lambda
#   # -0.00901  0.11646
#
# })

