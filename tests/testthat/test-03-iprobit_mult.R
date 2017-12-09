context("Multinomial estimation methods")

test_that("VB-EM (Closed-form)", {

  mod <- iprior::kernL(y ~ ., gen_mixture(n = 6, m = 3, seed = 123),
                       est.psi = FALSE)
  theta0 <- matrix(rep(1, 6), ncol = 3)
  alpha0 <- rep(1, 3)

  # Different hyperparameters in each class ------------------------------------
  mod1 <- iprobit(mod, control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                      alpha0 = alpha0))
  mod2 <- iprobit_mult(mod, 20, silent = TRUE, alpha0 = alpha0, theta0 = theta0)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(6.502061e-01, 5.067932e-01, -6.667536e-05), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(0.998459716, 0.001371398, 0.760741027), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(1.3513342, 0.1209030, 0.6443124), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

  # Only different RKHS parameters in each class -------------------------------
  mod1 <- iprobit(mod, common.intercept = TRUE,
                  control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                 alpha0 = alpha0))
  mod2 <- iprobit_mult(mod, 20, silent = TRUE, alpha0 = alpha0, theta0 = theta0,
                       common.intercept = TRUE)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(1, 0.4390309548, -0.0001031719), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(1, 0.0009601449, 0.7745982109), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(1, 0.1150476, 0.6924896), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

  # Only different intercepts in each class ------------------------------------
  mod1 <- iprobit(mod, common.RKHS.scale = TRUE,
                  control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                 alpha0 = alpha0))
  mod2 <- iprobit_mult(mod, 20, silent = TRUE, alpha0 = alpha0, theta0 = theta0,
                       common.RKHS.scale = TRUE)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(0.61014, -0.005820041, 0.048168828), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(0.99696, -0.005820041, 0.048168828), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(1.39291, -0.005820041, 0.048168828), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

  # Same hyperparameters in each class -----------------------------------------
  mod1 <- iprobit(mod, common.intercept = TRUE, common.RKHS.scale = TRUE,
                  control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                 alpha0 = alpha0))
  mod2 <- iprobit_mult(mod, 20, silent = TRUE, alpha0 = alpha0, theta0 = theta0,
                       common.intercept = TRUE, common.RKHS.scale = TRUE)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(1, 0 , 0), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

})

test_that("VB-EM with Metropolis sampler", {

  mod <- iprior::kernL(y ~ ., gen_circle(n = 6, m = 3, seed = 123), kernel = "fbm",
                       est.hurst = TRUE, est.psi = FALSE)
  theta0 <- matrix(c(1, 1, 0, 0), ncol = 3, nrow = 4)
  alpha0 <- rep(1, 3)

  # Different hyperparameters in each class ------------------------------------
  mod1 <- iprobit(mod, control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                      alpha0 = alpha0, n.samp = 5,
                                      thin.samp = 1, seed = 123))
  mod2 <- iprobit_mult_metr(mod, 20, silent = TRUE, alpha0 = alpha0,
                            theta0 = theta0, n.samp = 5, thin.samp = 1,
                            seed = 123)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(0.9985691, 1.0150100, 1.1014431, 0.5292299, 0.6026026),
               tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(1.0055504, 0.6488293, 0.9076343, 0.5466233, 0.4670734),
               tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(0.9958806, 0.9745103, 0.9592181, 0.4971903, 0.4359529),
               tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

  # Only different RKHS parameters in each class -------------------------------
  mod1 <- iprobit(mod, common.intercept = TRUE,
                  control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                 alpha0 = alpha0, n.samp = 5, thin.samp = 1,
                                 seed = 123))
  mod2 <- iprobit_mult_metr(mod, 20, silent = TRUE, alpha0 = alpha0,
                            theta0 = theta0, n.samp = 5, thin.samp = 1,
                            seed = 123, common.intercept = TRUE)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(1, 1.0150100, 1.1014431, 0.5292299, 0.6026026), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(1, 0.6488293, 0.9076343, 0.5466233, 0.4670734), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(1, 0.9745103, 0.9592181, 0.4971903, 0.4359529), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

  # Only different intercepts in each class ------------------------------------
  mod1 <- iprobit(mod, common.RKHS.scale = TRUE,
                  control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                 alpha0 = alpha0, n.samp = 5, thin.samp = 1,
                                 seed = 123))
  mod2 <- iprobit_mult_metr(mod, 20, silent = TRUE, alpha0 = alpha0,
                            theta0 = theta0, n.samp = 5, thin.samp = 1,
                            seed = 123, common.RKHS.scale = TRUE)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(0.9984867, 1.0390585, 1.0274495, 0.4885230, 0.5595608),
               tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(1.0051488, 1.0390585, 1.0274495, 0.4885230, 0.5595608),
               tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(0.9963646, 1.0390585, 1.0274495, 0.4885230, 0.5595608),
               tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

  # Same hyperparameters in each class -----------------------------------------
  mod1 <- iprobit(mod, common.intercept = TRUE, common.RKHS.scale = TRUE,
                  control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                 alpha0 = alpha0, n.samp = 5, thin.samp = 1,
                                 seed = 123))
  mod2 <- iprobit_mult_metr(mod, 20, silent = TRUE, alpha0 = alpha0,
                            theta0 = theta0, n.samp = 5, thin.samp = 1,
                            seed = 123, common.RKHS.scale = TRUE,
                            common.intercept = TRUE)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(1, 1.0390585, 1.0274495, 0.4885230, 0.5595608), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

})

context("Multinomial models")

# Generate data and model fit --------------------------------------------------
set.seed(123)
dat <- gen_circle(n = 15, m = 3)
dat.test <- gen_circle(n = 6, m = 3)

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
  expect_s3_class(mod, "iprobitMod_mult")
  expect_s3_class(modf, "iprobitMod")
  expect_s3_class(modf, "iprobitMod_mult")

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
