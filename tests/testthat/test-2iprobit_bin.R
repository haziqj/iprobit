context("Binary models")

test_that("Fitting binary models", {

	dat <- gen_mixture(n = 10)
	mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
	expect_s3_class(mod, "iprobitMod")
	expect_s3_class(mod, "iprobitMod_bin")

})

test_that("Print and summary", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  expect_that(print(mod), prints_text("Lower bound value ="))
  mod.summary <- summary(mod)
  expect_s3_class(mod.summary, "iprobitSummary")

})

test_that("Fitted and predict", {

  dat <- gen_mixture(n = 10)
  dat.test <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  expect_that(fitted(mod), is_a("list"))
  mod.predict <- predict(mod, newdata = dat.test$X)
  expect_s3_class(mod.predict, "iprobitPredict")
  expect_that(print(mod.predict), prints_text("Test data not provided."))
  mod.predict <- predict(mod, newdata = dat.test$X, y = dat.test$y)
  expect_that(print(mod.predict), prints_text("Test error rate"))

})

test_that("Convergence", {

  set.seed(123)
  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, control = list(maxit = 500, silent = TRUE))
  expect_equal(mod$lambda, 0.11646, tolerance = 1e-4)

  # Single lambda
  # mod <- iprobit_bin(dat$y, dat$X, silent = TRUE, maxit = 200)
  # > mod
  # Lower bound value =  -4.89696
  # Iterations =  106
  #
  # alpha   lambda
  # -0.00901  0.11646

})
