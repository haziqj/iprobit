data(iris)
# iris <- iris[c(1:3, 51:53, 101:103), ]
y <- iris$Species
y.lev <- levels(y)
X <- iris[, -5]

mod <- iprobit_mult(y, X)
mod$lower.bound[length(mod$lower.bound)]
