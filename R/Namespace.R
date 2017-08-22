# iprobit package imports ------------------------------------------------------
#' @importFrom stats coef delete.response dnorm fitted integrate kernel
#'   model.extract model.frame pnorm predict qnorm rnorm logLik
#' @importFrom utils capture.output data head setTxtProgressBar txtProgressBar
#' @import iprior
#' @import ggplot2
NULL

# Hacky way to pass R CMD CHECK "no visible binding" note ----------------------
globalVariables(c(".hMatList", "BlockBstuff", "Class", "H2l", "Hl", "Hlam.mat",
                  "Hlam.matsq", "Hsql", "Hurst", "Iteration", "Pl", "Psql",
                  "Sl", "X.1", "X.2", "X1", "X2", "Y", "alpha0", "class2",
                  "common.RKHS.scale", "common.intercept", "i", "ind1", "ind2",
                  "intr", "intr.3plus", "ipriorKernel", "iprobit.env", "l",
                  "lambda", "lambda.sq", "lambda0", "m", "maxit", "mod",
                  "model", "n", "no.int", "one.lam", "stop.crit", "value",
                  "variable", "w", "w0", "x", "y.levels", "yname", "ystar"))
