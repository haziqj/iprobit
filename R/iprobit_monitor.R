iprobit_monitor <- function(formula, data, test.data = NULL, kernel = "Canonical",
                            one.lam = FALSE, parsm = TRUE, monitoring = 10) {
  for (i in seq_len(monitoring)) {
    mod <- iprobit(formula, data, kernel = kernel, one.lam = one.lam, parsm = parsm,
                   control = list(maxit = 1))
    p1 <- iplot_predict(mod)
    print(p1)
  }
}
