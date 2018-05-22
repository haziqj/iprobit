context("Binary estimation methods")

test_that("VB-EM (Closed-form)", {

  mod <- iprior::kernL(y ~ ., gen_mixture(n = 4, seed = 123), est.psi = FALSE)
  mod1 <- iprobit(mod, control = list(silent = TRUE, maxit = 20, theta0 = 1:2,
                                      alpha0 = 1))
  mod2 <- iprobit_bin(mod, 20, silent = TRUE, alpha0 = 1, theta0 = 1:2)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(-0.070530717, 0.147493418, 0.004843337), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

})

# test_that("VB-EM with Metropolis sampler", {
#
#   mod <- iprior::kernL(y ~ ., gen_circle(n = 4, seed = 123), kernel = "se",
#                        est.lengthscale = TRUE, est.psi = FALSE)
#   mod1 <- iprobit(mod, control = list(silent = TRUE, maxit = 2, theta0 = 1:4,
#                                       alpha0 = 1, n.samp = 10, thin.samp = 1,
#                                       seed = 123))
#   mod2 <- iprobit_bin_metr(mod, 2, silent = TRUE, alpha0 = 1, theta0 = 1:4,
#                            n.samp = 10, thin.samp = 1, seed = 123)
#   expect_equal(as.numeric(coef(mod1)[, 1]),
#                c(0.1395819, 1.0103551, 2.0785975, 18.2429103, 83.1534269),
#                tol = 1e-5)
#   expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)
#
# })

context("Binary models")

# Generate data and model fit --------------------------------------------------
dat <- gen_mixture(n = 10)
dat.test <- gen_mixture(n = 5)

# Regular
mod <- iprobit(dat$y, dat$X, control = list(maxit = 5, silent = TRUE))
# Formula (with training samples)
modf <- iprobit(y ~ ., dat, one.lam = TRUE, train.samp = c(1, 5),
                control = list(maxit = 5, silent = TRUE))
# Metropolis
modmetr <- iprobit(y ~ ., dat, kernel = "fbm", est.hurst = TRUE,
                   control = list(maxit = 2, silent = TRUE))

# Predict without quantiles ----------------------------------------------------
mod.predict <- predict(mod, newdata = list(dat.test$X))
modf.predict <- predict(modf, newdata = dat.test)
modmetr.predict <- predict(modmetr, newdata = dat.test)

# Predict with quantiles -------------------------------------------------------
mod.predict.q <- predict(mod, list(dat.test$X), quantiles = TRUE, n.samp = 3)
modf.predict.q <- predict(modf, dat.test, quantiles = TRUE, n.samp = 3)
modmetr.predict.q <- predict(modmetr, dat.test,  quantiles = TRUE, n.samp = 3)

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
  expect_s3_class(modmetr.predict, "iprobitPredict")

  expect_s3_class(mod.predict.q, "iprobitPredict_quant")
  expect_s3_class(modf.predict.q, "iprobitPredict_quant")
  expect_s3_class(modmetr.predict.q, "iprobitPredict_quant")

  expect_that(print(mod.predict), prints_text("Test data not provided."))

})

test_that("Predict (with test error rate)", {

  mod.predict <- predict(mod, newdata = list(dat.test$X), y = dat.test$y)

  expect_that(print(mod.predict), prints_text("Test error"))
  expect_that(print(modf.predict), prints_text("Test error"))
  expect_that(print(modmetr.predict), prints_text("Test error"))

})

# test_that("Convergence", {
#
#   set.seed(123)
#   dat <- gen_mixture(n = 10)
#   mod <- iprobit(dat$y, dat$X, control = list(maxit = 500, silent = TRUE))
#   modf <- iprobit(y ~ ., dat, one.lam = TRUE,
#                   control = list(maxit = 500, silent = TRUE))
#   expect_equal(as.numeric(get_lambda(mod)), 0.11646, tolerance = 1e-3)
#   expect_equal(as.numeric(get_lambda(mod)), 0.11646, tolerance = 1e-3)
#
#   # > summary(mod)
#   # Call:
#   # iprobit(y = dat$y, X1 = dat$X, control = list(maxit = 500, silent = TRUE))
#   #
#   # Classes:
#   # RKHS used:
#   # Linear (X1)
#   #
#   # Hyperparameters:
#   #           Mean   S.D.    2.5%  97.5%
#   # alpha  -0.0090 0.3162 -0.6288 0.6108
#   # lambda  0.1165 0.0201  0.0770 0.1559
#   # ---
#   #
#   # Closed-form VB-EM algorithm. Iterations: 71/500
#   # Converged to within 1e-05 tolerance. Time taken: 0.08387208 secs
#   # Variational lower bound: -4.897003
#   # Training error: 0%. Brier score: 0.004980669
#
# })
