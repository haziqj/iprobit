# This is wrong
# EprodPhiZ_1 <- function(mu, sigma) {
#   Sigma <- outer(sigma, sigma, FUN = "*") + diag(1, length(mu))
#   L.inv <- solve(t(chol(Sigma)))
#   e <- L.inv %*% mu
#   prod(pnorm(e))
# }
#
# EprodPhiZ_2 <- function(mu, sigma) {
#   Sigma <- outer(sigma, sigma, FUN = "*") + diag(1, length(mu))
#   pmvnorm(upper = mu, sigma = Sigma)
# }

EprodPhiZ <- function(mu, sigma = rep(1, length(mu)), log = FALSE) {
  res <- integrate(
    function(z) {
      tmp <- 0
      for (i in seq_len(length(mu)))
        tmp <- tmp + pnorm(z + mu[i], log.p = TRUE)
      exp(tmp) * dnorm(z)
    },
    lower = -Inf, upper = Inf
  )$value
  ifelse(isTRUE(log), log(res), res)
}

# EprodPhiZ <- function(mu, sigma) {
#   res1 <- EprodPhiZ_1(mu = mu, sigma = sigma)
#   res2 <- EprodPhiZ_2(mu = mu, sigma = sigma)
#   res3 <- EprodPhiZ_3(mu = mu, sigma = sigma)
#
#   list(res1, res2, res3)
# }
#
# k <- 20
# microbenchmark::microbenchmark(
#   method1 = EprodPhiZ_1(rnorm(k), rep(1, k)),
#   method2 = EprodPhiZ_2(rnorm(k), rep(1, k)),
#   method3 = EprodPhiZ_3(rnorm(k), rep(1, k))
# )
#
#
# EprodPhiZ(rnorm(4), rep(1, 4))

