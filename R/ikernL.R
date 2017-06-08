#' @export
ikernL <- function(Xl, newdata = NULL, kernel = c("Canonical", "FBM,0.5", "Pearson"),
                   interactions = NULL) {
  Hurst <- splitHurst(kernel)  # get the Hurst coefficient
  Hurst <- ifelse(is.na(Hurst), 0.5, Hurst)
  kernel <- splitKernel(kernel)  # get the kernel
  kernel <- match.arg(kernel, c("Canonical", "FBM", "Pearson"))
  if (kernel == "FBM")
    kernelFn <- function(x, y = NULL) iprior::fnH3(x = x, y = y, gamma = Hurst)
  else
    kernelFn <- function(x, y = NULL) iprior::fnH2(x = x, y = y)
  p <- length(Xl)
  Hl <- NULL

  for (i in 1:p) {
    if (is.factor(Xl[[i]]))
      Hl[[i]] <- iprior::fnH1(x = Xl[[i]], y = newdata[[i]])
    else
     Hl[[i]] <- kernelFn(Xl[[i]], y = newdata[[i]])
  }
  Hl
}

# # Test iprior kernel loader on categorical data
#
# dat <- gen_mixture(n = 2)
# df <- data.frame(dat$X, y = dat$y)
# mod <- kernL(y ~ . ^ 2, df)
#
# myfun <- function() {
#   this.env <- environment()
#   environment(.lambdaExpand) <- this.env
#   list2env(mod, this.env)
#   list2env(BlockBstuff, this.env)
#   environment(BlockB)  <- this.env
#   list2env(model, this.env)
#   lambda <- seq_len(p)
#   .lambdaExpand(env = this.env)
#   BlockB(1)
#   BlockB(2)
#   print(lambda)
#   print(Hl)
#   print(Pl)
#   print(Psql)
#   print(Sl)
# }
# myfun()




