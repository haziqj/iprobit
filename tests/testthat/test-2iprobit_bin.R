context("Binary models")

dat <- gen_mixture(n = 10)
dat.test <- gen_mixture(n = 5)

mod <- iprobit(dat$y, dat$X, control = list(maxit = 5, silent = TRUE))
modf <- iprobit(y ~ ., dat, one.lam = TRUE, control = list(maxit = 5,
                                                           silent = TRUE))

dat.test <- gen_mixture(n = 5)
mod.predict <- predict(mod, newdata = list(dat.test$X))
modf.predict <- predict(modf, newdata = dat.test)

mod.predict.q <- predict(mod, newdata = list(dat.test$X), quantiles = TRUE,
                         n.samp = 3)
modf.predict.q <- predict(modf, newdata = dat.test, quantiles = TRUE,
                          n.samp = 3)


test_that("Fitting binary models", {

	expect_s3_class(mod, "iprobitMod")
	expect_s3_class(mod, "iprobitMod_bin")
	expect_s3_class(modf, "iprobitMod")
	expect_s3_class(modf, "iprobitMod_bin")

})

test_that("Print", {

  expect_that(print(mod), prints_text("Training error rate"))
  expect_that(print(modf), prints_text("Training error rate"))

})

test_that("Summary", {

  expect_s3_class(summary(mod), "iprobitMod_summary")
  expect_s3_class(summary(modf), "iprobitMod_summary")

})

test_that("Fitted", {

  expect_that(fitted(mod), is_a("iprobitPredict"))
  expect_that(fitted(modf), is_a("iprobitPredict"))

})

test_that("Fitted quantiles", {

  expect_that(fitted(mod, TRUE, n.samp = 3), is_a("iprobitPredict_quant"))
  expect_that(fitted(modf, TRUE, n.samp = 3), is_a("iprobitPredict_quant"))

})

test_that("Predict (without test error rate)", {

  expect_s3_class(mod.predict, "iprobitPredict")
  expect_s3_class(modf.predict, "iprobitPredict")
  expect_s3_class(mod.predict.q, "iprobitPredict_quant")
  expect_s3_class(modf.predict.q, "iprobitPredict_quant")
  expect_that(print(mod.predict), prints_text("Test data not provided."))

})

test_that("Predict (with test error rate)", {

  mod.predict <- predict(mod, newdata = list(dat.test$X), y = dat.test$y)

  expect_that(print(mod.predict), prints_text("Test error"))
  expect_that(print(modf.predict), prints_text("Test error"))

})

test_that("Convergence", {

  set.seed(123)
  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, control = list(maxit = 500, silent = TRUE))
  modf <- iprobit(y ~ ., dat, one.lam = TRUE,
                  control = list(maxit = 500, silent = TRUE))
  expect_equal(as.numeric(get_lambda(mod)), 0.11646, tolerance = 1e-3)
  expect_equal(as.numeric(get_lambda(mod)), 0.11646, tolerance = 1e-3)

  # > summary(mod)
  # Call:
  # iprobit(y = dat$y, X1 = dat$X, control = list(maxit = 500, silent = TRUE))
  #
  # Classes:
  # RKHS used:
  # Linear (X1)
  #
  # Hyperparameters:
  #           Mean   S.D.    2.5%  97.5%
  # alpha  -0.0090 0.3162 -0.6288 0.6108
  # lambda  0.1165 0.0201  0.0770 0.1559
  # ---
  #
  # Closed-form VB-EM algorithm. Iterations: 71/500
  # Converged to within 1e-05 tolerance. Time taken: 0.08387208 secs
  # Variational lower bound: -4.897003
  # Training error: 0%. Brier score: 0.004980669

})
