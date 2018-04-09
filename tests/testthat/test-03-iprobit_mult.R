context("Multinomial estimation methods")

test_that("VB-EM (Closed-form)", {

  mod <- iprior::kernL(y ~ ., gen_mixture(n = 6, m = 3, seed = 123),
                       est.psi = FALSE)
  theta0 <- matrix(rep(1, 6), ncol = 3)
  alpha0 <- rep(1, 3)

  mod1 <- iprobit(mod, control = list(silent = TRUE, maxit = 20,
                                      theta0 = theta0, alpha0 = alpha0))
  mod2 <- iprobit_mult(mod, 20, silent = TRUE, alpha0 = alpha0, theta0 = theta0)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(-0.389863967, -0.005820041, 0.048168828), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(-0.003042599, -0.005820041, 0.048168828), tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(0.392906565, -0.005820041, 0.048168828), tol = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tol = 1e-5)

})

test_that("VB-EM with Metropolis sampler", {

  mod <- iprior::kernL(y ~ ., gen_circle(n = 6, m = 3, seed = 123), kernel = "fbm",
                       est.hurst = TRUE, est.psi = FALSE)
  theta0 <- matrix(c(1, 1, 0, 0), ncol = 3, nrow = 4)
  alpha0 <- rep(1, 3)

  # Only different intercepts in each class ------------------------------------
  mod1 <- iprobit(mod,
                  control = list(silent = TRUE, maxit = 20, theta0 = theta0,
                                 alpha0 = alpha0, n.samp = 5, thin.samp = 1,
                                 seed = 123))
  mod2 <- iprobit_mult_metr(mod, 20, silent = TRUE, alpha0 = alpha0,
                            theta0 = theta0, n.samp = 5, thin.samp = 1,
                            seed = 123)
  expect_equal(as.numeric(coef(mod1)[, 1]),
               c(-0.001513346, 1.0390585, 1.0274495, 0.4885230, 0.5595608),
               tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 2]),
               c(0.005148759, 1.0390585, 1.0274495, 0.4885230, 0.5595608),
               tol = 1e-5)
  expect_equal(as.numeric(coef(mod1)[, 3]),
               c(-0.003635412, 1.0390585, 1.0274495, 0.4885230, 0.5595608),
               tol = 1e-5)
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
