R/iprobit: Binary and multinomial probit regression using I-priors
================

This is an `R` package which extends I-prior regression to unordered categorical responses via a probit link function. This allows the user to fit models for classification or inference using fitted probabilities. Estimation is performed using a variational EM algorithm.

Binary classification (toy example)
-----------------------------------

#### Model fitting

``` r
dat <- gen_spiral(n = 300)  # generate binary toy example data set
mod <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "FBM")
## ========================================================================
## Converged after 96 iterations.
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
## lambda 5.6704 0.2322  5.2153 6.1254
## 
## Converged to within 1e-05 tolerance. No. of iterations: 96
## Variational lower bound: -140.7216
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
## Lower bound value =  -258.8624 
## Iterations =  10 
## 
##            Class = 1 Class = 2 Class = 3 Class = 4
## alpha       -0.31319  -0.72351  -0.51619  -0.65401
## lambda[1,]   1.62060   0.00000   0.04799   0.00000
## lambda[2,]   0.00000   0.30576   0.00000   0.65504
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
## Test error rate: 7 %

mod.pred$y
##   [1] 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2
##  [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
##  [71] 3 3 3 3 4 4 4 4 3 1 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3
## Levels: 1 2 3 4

head(mod.pred$prob)
##        1      2      3      4
## 1 0.7933 0.2014 0.0015 0.0038
## 2 0.8121 0.0208 0.0038 0.1633
## 3 0.9095 0.0768 0.0007 0.0130
## 4 0.7169 0.0815 0.0287 0.1729
## 5 0.4643 0.2842 0.1476 0.1039
## 6 0.9658 0.0252 0.0000 0.0090
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
## Variational lower bound: -50.48093
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
