context("Binary models")

test_that("Fitting binary models", {

	dat <- gen_mixture(n = 10)
	mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
	modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE, control = list(maxit = 5))
	expect_s3_class(mod, "iprobitMod")
	expect_s3_class(mod, "iprobitMod_bin")
	expect_s3_class(modf, "iprobitMod")
	expect_s3_class(modf, "iprobitMod_bin")

})

test_that("Print", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE, control = list(maxit = 5))
  expect_that(print(mod), prints_text("Lower bound value ="))
  expect_that(print(modf), prints_text("Lower bound value ="))
  mod.summary <- summary(mod)
  expect_s3_class(mod.summary, "iprobitSummary")

})

test_that("Summary", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE, control = list(maxit = 5))
  mod.summary <- summary(mod)
  modf.summary <- summary(modf)
  expect_s3_class(mod.summary, "iprobitSummary")
  expect_s3_class(modf.summary, "iprobitSummary")

})

test_that("Fitted", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE, control = list(maxit = 5))
  expect_that(fitted(mod), is_a("list"))
  expect_that(fitted(modf), is_a("list"))

})

test_that("Predict (without test error rate)", {

  dat <- gen_mixture(n = 10)
  dat.test <- gen_mixture(n = 5)
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

  dat <- gen_mixture(n = 10)
  dat.test <- gen_mixture(n = 5)
  mod <- iprobit(dat$y, dat$X, silent = TRUE, control = list(maxit = 5))
  modf <- iprobit(y ~ ., dat, silent = TRUE, one.lam = TRUE,
                  control = list(maxit = 5))

  mod.predict <- predict(mod, newdata = list(dat.test$X), y = dat.test$y)
  modf.predict <- predict(modf, newdata = dat.test)

  expect_that(print(mod.predict), prints_text("Test error rate"))
  expect_that(print(modf.predict), prints_text("Test error rate"))

})

test_that("Convergence", {

  set.seed(123)
  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, control = list(maxit = 500, silent = TRUE))
  modf <- iprobit(y ~ ., dat, one.lam = TRUE, silent = TRUE)
  expect_equal(mod$lambda, 0.11646, tolerance = 1e-3)
  expect_equal(modf$lambda, 0.11646, tolerance = 1e-3)

  # Single lambda
  # mod <- iprobit_bin(dat$y, dat$X, silent = TRUE, maxit = 200)
  # > mod
  # Lower bound value =  -4.89696
  # Iterations =  106
  #
  # alpha   lambda
  # -0.00901  0.11646

})
