context("Kernels and interactions")

# Binary probit models ---------------------------------------------------------

test_that("FBM kernel (binary)", {

  dat <- gen_mixture(n = 10)

  # FBM = 0.5
  mod <- iprobit(y ~ ., dat, kernel = "FBM", silent = TRUE)
  expect_equivalent(get_Hurst(mod), c(0.5, 0.5))
  mod <- iprobit(dat$y, dat$X, kernel = "FBM", silent = TRUE)
  expect_equivalent(get_Hurst(mod), 0.5)

  # FBM = 0.1, 0.9
  mod <- iprobit(y ~ ., dat, kernel = c("FBM,0.1", "FBM,0.9"), silent = TRUE)
  expect_equivalent(get_Hurst(mod), c(0.1, 0.9))

})

test_that("Pearson kernel (binary)", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$y, silent = TRUE)
  expect_equivalent(get_kernel(mod), "Pearson")

})

test_that("Mixed kernel (binary)", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit(y ~ ., dat, kernel = c("Canonical", "FBM"), silent = TRUE)
  expect_equivalent(get_kernel(mod), c("Canonical", "FBM,0.5"))

})

test_that("Single lambda (binary)", {

  dat <- gen_mixture(n = 10)
  mod <- iprobit(dat$y, dat$X, silent = TRUE,
                 control = list(alpha0 = 1, lambda0 = rep(1, 2)))
  modf <- iprobit(y ~ ., dat, one.lam = TRUE, silent = TRUE,
                  control = list(alpha0 = 1, lambda0 = rep(1, 2)))
  expect_equivalent(coef(mod), coef(modf))

})

test_that("Parsimonious interactions (binary)", {

  dat <- gen_mixture(n = 100)
  mod <- iprobit(dat$y, dat$X[, 1], dat$X[, 2], silent = TRUE,
                 interactions = "1:2", parsm = TRUE,
                 control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2)))
  modf <- iprobit(y ~ . ^ 2, dat, silent = TRUE, parsm = TRUE,
                  control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2)))
  expect_equivalent(coef(mod), coef(modf))

})

test_that("Non-parsimonious interactions (binary)", {

  dat <- gen_mixture(n = 100)
  mod <- iprobit(dat$y, dat$X[, 1], dat$X[, 2], silent = TRUE,
                 interactions = "1:2", parsm = FALSE,
                 control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 3)))
  modf <- iprobit(y ~ . ^ 2, dat, silent = TRUE, parsm = FALSE,
                  control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 3)))
  expect_equivalent(coef(mod), coef(modf))

})

test_that("Squared terms (binary)", {

  dat <- gen_mixture(n = 100)
  mod <- iprobit(dat$y, dat$X[, 1], dat$X[, 1] ^ 2, silent = TRUE,
                 control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2)))
  modf <- iprobit(y ~ X1 + I(X1 ^ 2), dat, silent = TRUE,
                  control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2)))
  expect_equivalent(coef(mod), coef(modf))

})

# Multinomial probit models ----------------------------------------------------

test_that("FBM kernel (multinomial)", {

  dat <- gen_mixture(n = 3, m = 3)

  # FBM = 0.5
  mod <- iprobit(y ~ ., dat, kernel = "FBM", silent = TRUE)
  expect_equivalent(get_Hurst(mod), c(0.5, 0.5))
  mod <- iprobit(dat$y, dat$X, kernel = "FBM", silent = TRUE)
  expect_equivalent(get_Hurst(mod), 0.5)

  # FBM = 0.1, 0.9
  mod <- iprobit(y ~ ., dat, kernel = c("FBM,0.1", "FBM,0.9"), silent = TRUE)
  expect_equivalent(get_Hurst(mod), c(0.1, 0.9))

})

test_that("Pearson kernel (multinomial)", {

  dat <- gen_mixture(n = 3, m = 3)
  mod <- iprobit(dat$y, dat$y, silent = TRUE)
  expect_equivalent(get_kernel(mod), "Pearson")

})

test_that("Mixed kernel (multinomial)", {

  dat <- gen_mixture(n = 3, m = 3)
  mod <- iprobit(y ~ ., dat, kernel = c("Canonical", "FBM"), silent = TRUE)
  expect_equivalent(get_kernel(mod), c("Canonical", "FBM,0.5"))

})

test_that("Single lambda (multinomial)", {

  dat <- gen_mixture(n = 3, m = 3)
  mod <- iprobit(dat$y, dat$X, silent = TRUE,
                 control = list(alpha0 = 1, lambda0 = rep(1, 2 * 3)))
  modf <- iprobit(y ~ ., dat, one.lam = TRUE, silent = TRUE,
                  control = list(alpha0 = 1, lambda0 = rep(1, 2 * 3)))
  expect_equivalent(coef(mod), coef(modf))

})

test_that("Parsimonious interactions (multinomial)", {

  dat <- gen_mixture(n = 100, m = 3)
  mod <- iprobit(dat$y, dat$X[, 1], dat$X[, 2], silent = TRUE,
                 interactions = "1:2", parsm = TRUE,
                 control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2 * 3)))
  modf <- iprobit(y ~ . ^ 2, dat, silent = TRUE, parsm = TRUE,
                  control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2 * 3)))
  expect_equivalent(coef(mod), coef(modf))

})

test_that("Non-parsimonious interactions (multinomial)", {

  dat <- gen_mixture(n = 100, m = 3)
  mod <- iprobit(dat$y, dat$X[, 1], dat$X[, 2], silent = TRUE,
                 interactions = "1:2", parsm = FALSE,
                 control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 3 * 3)))
  modf <- iprobit(y ~ . ^ 2, dat, silent = TRUE, parsm = FALSE,
                  control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 3 * 3)))
  expect_equivalent(coef(mod), coef(modf))

})

test_that("Squared terms (multinomial)", {

  dat <- gen_mixture(n = 100, m = 3)
  mod <- iprobit(dat$y, dat$X[, 1], dat$X[, 1] ^ 2, silent = TRUE,
                 control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2)))
  modf <- iprobit(y ~ X1 + I(X1 ^ 2), dat, silent = TRUE,
                  control = list(maxit = 3, alpha0 = 1, lambda0 = rep(1, 2)))
  expect_equivalent(coef(mod), coef(modf))

})


