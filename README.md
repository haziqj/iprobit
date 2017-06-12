R/iprobit: Binary and multinomial probit regression using I-priors
================

[![Build Status](https://travis-ci.org/haziqjamil/iprobit.svg?branch=master)](https://travis-ci.org/haziqjamil/iprobit) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/haziqjamil/iprobit?branch=master&svg=true)](https://ci.appveyor.com/project/haziqjamil/iprobit) [![Coverage Status](https://img.shields.io/codecov/c/github/haziqjamil/iprobit/master.svg)](https://codecov.io/gh/haziqjamil/iprobit)

This is an `R` package which extends I-prior regression to unordered categorical responses via a probit link function. This allows the user to fit models for classification or inference using fitted probabilities. Estimation is performed using a variational EM algorithm.

Binary classification (toy example)
-----------------------------------

#### Model fitting

``` r
dat <- gen_spiral(n = 300)  # generate binary toy example data set
mod <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "FBM")
## ===========================================================================
## Convergence criterion not met.
```

#### Model summary

``` r
summary(mod)
## 
## Call:
## iprobit(formula = y ~ X1 + X2, data = dat)
## 
## Classes: 1, 2 
## 
## RKHS used:
## Fractional Brownian motion with Hurst coef. 0.5 (X1 + X2) 
## 
## Parameter estimates:
##          Mean   S.D.    2.5%  97.5%
## alpha  0.0000 0.0577 -0.1132 0.1132
## lambda 5.6714 0.2321  5.2164 6.1265
## 
## Convergence criterion not met. No. of iterations: 100
## Variational lower bound: -140.7163
```

#### Boundary plot for two-dimensional covariates

``` r
iplot_predict(mod)
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

Multiclass classification (toy example)
---------------------------------------

#### Model fit report and parameter estimates

``` r
dat <- gen_mixture(n = 500, m = 4, sd = 1.5)  # generate 4-class toy example data set
(mod <- iprobit(y ~ X1 + X2, dat, control = list(maxit = 10)))
## ===========================================================================
## Convergence criterion not met.
## Lower bound value =  -257.5504 
## Iterations =  10 
## 
##            Class = 1 Class = 2 Class = 3 Class = 4
## alpha        0.01109  -0.01772  -0.05802  -0.15256
## lambda[1,]   0.65506   0.00000   1.11217   0.00000
## lambda[2,]   0.00000   1.37593   0.00000   0.71118
```

#### Boundary plot for two-dimensional covariates

``` r
iplot_predict(mod)
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

#### Obtain out-of-sample test error rates, predicted classes and probabilities

``` r
dat.test <- gen_mixture(n = 100, m = 4, sd = 1.5)
(mod.pred <- predict(mod, newdata = dat.test))
## Test error rate: 8 %
## 
## Predicted classes:
##   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2
##  [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1 3 3 3
##  [71] 3 3 3 3 3 3 3 4 4 4 4 1 3 4 4 4 3 4 4 4 4 4 4 4 4 4 3 4 4 4
## Levels: 1 2 3 4
## 
## Predicted probabilities:
##         1      2      3      4
## 1  0.9594 0.0240 0.0000 0.0166
## 2  0.5279 0.4386 0.0156 0.0179
## 3  0.9557 0.0167 0.0000 0.0276
## 4  0.8895 0.0761 0.0003 0.0341
## 5  0.4983 0.0250 0.0201 0.4566
## 6  0.6884 0.0916 0.0140 0.2060
## 7  0.9827 0.0150 0.0000 0.0023
## 8  0.5008 0.4945 0.0026 0.0021
## 9  0.8401 0.0900 0.0015 0.0684
## 10 0.7881 0.0529 0.0028 0.1562
## ...

mod.pred$y
##   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2
##  [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1 3 3 3
##  [71] 3 3 3 3 3 3 3 4 4 4 4 1 3 4 4 4 3 4 4 4 4 4 4 4 4 4 3 4 4 4
## Levels: 1 2 3 4

head(mod.pred$prob)
##        1      2      3      4
## 1 0.9594 0.0240 0.0000 0.0166
## 2 0.5279 0.4386 0.0156 0.0179
## 3 0.9557 0.0167 0.0000 0.0276
## 4 0.8895 0.0761 0.0003 0.0341
## 5 0.4983 0.0250 0.0201 0.4566
## 6 0.6884 0.0916 0.0140 0.2060
```

Fisher's Iris data set
----------------------

#### Model fitting (common RKHS scale across classes for each covariate)

``` r
mod <- iprobit(Species ~ ., iris, kernel = "FBM", one.lam = TRUE,
               control = list(alpha0 = 1, lambda0 = 1, 
                              stop.crit = 1e-1,
                              common.RKHS.scale = TRUE, 
                              common.intercept = FALSE))
## ====================
## Converged after 27 iterations.

summary(mod)
## 
## Call:
## iprobit(formula = Species ~ ., data = iris)
## 
## Classes: setosa, versicolor, virginica 
## 
## RKHS used:
## Fractional Brownian motion with Hurst coef. 0.5 (Sepal.Length + ... + Petal.Width) 
## 
## Parameter estimates:
##            Mean   S.D.   2.5%  97.5%
## alpha[1] 0.8835 0.0816 0.7234 1.0435
## alpha[2] 1.0572 0.0816 0.8971 1.2172
## alpha[3] 1.0594 0.0816 0.8994 1.2194
## lambda   0.3474 0.0116 0.3246 0.3703
## 
## Converged to within 0.1 tolerance. No. of iterations: 27
## Variational lower bound: -50.38115
```

#### Monitor convergence

``` r
iplot_lb(mod)
```

![](README_files/figure-markdown_github/unnamed-chunk-8-1.png)

#### Plot of fitted probabilities

``` r
iplot_fitted(mod)
```

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

------------------------------------------------------------------------

Copyright (C) 2017 [Haziq Jamil](http://haziqj.ml).
