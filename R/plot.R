#' @import ggplot2
#' @export
plot.iprobitMod <- function(object, niter.plot = NULL, levels = NULL, ...) {
  iplot_fitted(object)
}

#' @export
iplot_fitted <- function(object) {
  list2env(object, environment())
  list2env(ipriorKernel, environment())
  list2env(model, environment())

  probs <- fitted(object)$prob
  df.plot <- data.frame(probs, i = 1:n)
  colnames(df.plot) <- c(y.levels, "i")
  df.plot <- reshape2::melt(df.plot, id.vars = "i")

  ggplot(df.plot, aes(x = i, y = value)) +
    geom_area(aes(col = variable, fill = variable), position = "stack",
              alpha = 0.95) +
    labs(col = "Class", fill = "Class", x = "Index", y = "Fitted probabilities") +
    coord_cartesian(expand = FALSE) +
    theme_bw()
}

# iplot_fitted_mult <- function(x, levels = NULL) {
#   list2env(x, envir = environment())
#   n <- length(y)
#   y.lev <- levels(y)
#   m <- length(y.lev)
#   nm <- n * m
#   p <- ncol(X)
#
#   probs <- fitted(x)$probs
#   df.plot <- data.frame(probs, i = 1:n)
#   df.plot <- reshape2::melt(df.plot, id.vars = "i")
#   ggplot(df.plot, aes(x = i, y = value)) +
#     geom_area(aes(col = variable, fill = variable), position = "stack") +
#     coord_cartesian(ylim = c(0, 1)) +
#     labs(col = "Class", fill = "Class", x = "Index", y = "Fitted probabilities") +
#     theme_bw()
# }

#' @export
iplot_lb <- function(x, niter.plot = NULL, lab.pos = c("up", "down")) {
  lab.pos <- match.arg(lab.pos, c("up", "down"))
  if (lab.pos == "up") lab.pos <- -0.5
  else lab.pos <- 1.5

  lb.original <- x$lower.bound[!is.na(x$lower.bound)]
  if (is.null(niter.plot)) niter.plot <- c(1, length(lb.original))
  else if (length(niter.plot) == 1) niter.plot <- c(1, niter.plot)
  niter.plot <- niter.plot[1]:niter.plot[2]
  lb <- lb.original[niter.plot]
  plot.df <- data.frame(Iteration = niter.plot, lb = lb)
  time.per.iter <- mod$time$time / (x$niter - 1)
  if (time.per.iter < 0.001) time.per.iter <- 0.001

  ggplot(plot.df, aes(x = Iteration, y = lb, label = max(lb))) +
    geom_line(col = "grey60") +
    geom_point() +
    geom_hline(yintercept = max(lb.original), linetype = 2, col = "red") +
    scale_x_continuous(
      sec.axis = sec_axis(~ . * time.per.iter, name = "Time (seconds)"),
      breaks = scales::pretty_breaks()
    ) +
    # directlabels::geom_dl(method = "bottom.pieces") +
    geom_text(aes(label = ifelse(lb == max(lb), round(max(lb), 2), "")),
              vjust = 1.5, size = 3.7) +
    annotate("text", col = "red", x = niter.plot[1], y = max(lb.original),
             vjust = lab.pos, label = round(max(lb.original), 2), size = 3.7) +
    labs(y = "Variational lower bound") +
    theme_bw()
}

#' @export
iplot_prob <- function(x, covariate = 1, levels = NULL) {
  if (!is.numeric(covariate)) {
    covariate <- grep(covariate, colnames(mod$X))
  }
  classes <- as.factor(x$y)
  if (!is.null(levels)) levels(classes) <- levels
  else levels(classes) <- x$y.levels
  mean.X <- apply(x$X, 2, mean)
  X.new <- x$X
  X.new[, -covariate] <- mean.X[-covariate]
  p.hat <- predict(x, X.new)$prob
  p.hat.lower <- predict(x, X.new, "lower")$prob
  p.hat.upper <- predict(x, X.new, "upper")$prob
  plot.df <- data.frame(x = x$X[, covariate], y = x$y, prob = p.hat,
                        low = p.hat.lower, high = p.hat.upper,
                        Class = classes)
  x.axis.lab <- colnames(x$X)[covariate]

  ggplot(plot.df, aes(x = x, y = y)) +
    geom_line(aes(x = x, y = prob)) +
    geom_point(aes(col = Class)) +
    geom_ribbon(aes(x = x, ymin = low, ymax = high), fill = "grey70", alpha = 0.5) +
    labs(x = x.axis.lab, y = "Probability") +
    theme_bw()
}

create_x_range <- function(X, X.var) {

}

#' @export
iplot_predict <- function(object, test.data = NULL, X.var = c(1, 2)) {
  X <- object$ipriorKernel$x
  class(X) <- NULL
  X <- as.data.frame(X)
  p <- ncol(X)
  xname <- object$ipriorKernel$model$xname[X.var]

  maxmin <- cbind(apply(X, 2, min), apply(X, 2, max))
  xx <- list(NULL)
  for (j in 1:2) {
    mm <- maxmin[X.var[j], ]
    xx[[j]] <- seq(from = mm[1] - 1, to = mm[2] + 1, length.out = 100)
  }
  mm <- maxmin[X.var, ]
  plot.df <- expand.grid(xx[[1]], xx[[2]])
  if (p > 2) {
    X.avg <- X[, -X.var]
    xx.avg <- apply(X.avg, 2, mean)
    plot.df <- cbind(plot.df, lapply(xx.avg, function(x) x))
  }

  classes <- factor(object$ipriorKernel$Y)
  levels(classes) <- object$ipriorKernel$y.levels
  points.df <- data.frame(X[, X.var], classes)

  if (!is.null(object$formula)) {
    # Fitted using formula
    xname <- colnames(plot.df)[1:2] <-
      attr(object$ipriorKernel$terms, "term.labels")[X.var]
    prob <- predict(object, newdata = plot.df)$prob
    if (!is.null(test.data)) {
      if (is.iprobitData(test.data)) test.data <- as.data.frame(test.data)
      points.df <- test.data
    }
  } else {
    prob <- predict(object, newdata = list(plot.df))$prob
    if (!is.null(test.data)) {
      if (is.iprobitData(test.data)) test.data <- as.data.frame(test.data)
      points.df <- test.data
    }
  }
  plot.df <- cbind(plot.df, prob)
  colnames(plot.df)[1:2] <- c("X1", "X2")
  colnames(plot.df)[-(1:p)] <- paste0("class", seq_along(object$y.levels))
  colnames(points.df) <- c("X1", "X2", "Class")

  if (is.iprobitMod_bin(object))
    p <- iplot_predict_bin(plot.df, points.df, mm[1, ], mm[2, ],
                           length(levels(classes)))
  if (is.iprobitMod_mult(object))
    p <- iplot_predict_mult(plot.df, points.df, mm[1, ],mm[2, ],
                            length(levels(classes)))

  yname <- ifelse(object$ipriorKernel$model$yname == "y", "Class",
                  object$ipriorKernel$model$yname)
  p <- p +
    labs(x = xname[1], y = xname[2], col = yname)
  p
}

iplot_predict_bin <- function(plot.df, points.df, x, y, m) {
  ggplot() +
    geom_raster(data = plot.df, aes(X1, X2, fill = class2), alpha = 0.5) +
    scale_fill_gradient(low = "#F8766D", high = "#00BFC4", limits = c(0, 1)) +
    geom_point(data = points.df, aes(X1, X2, col = Class)) +
    coord_cartesian(xlim = x, ylim = y) +
    guides(fill = FALSE) +
    theme_bw()
}

iplot_predict_mult <- function(plot.df, points.df, x, y, m) {
  fill.col <- iprior::ggColPal(m)
  alpha <- 0.6
  class.ind <- sort(seq(from = ncol(plot.df), by = -1, length = m))

  # decbound.df <- NULL
  # for (j in 1:m) {
  #   tmp <- plot.df[round(plot.df[, 2 + j], 1) == 0.5, 1:2]
  #   decbound.df[[j]] <- tmp[order(tmp[, 1]), ]
  # }

  # Add first layer ------------------------------------------------------------
  p <- ggplot() +
    geom_raster(data = plot.df, aes(X1, X2, fill = fill.col[1],
                                    alpha = alpha * plot.df[, class.ind[1]])) +
    scale_alpha_continuous(range = c(0, 0.5))
    # geom_line(data = decbound.df[[1]], aes(X1, X2))

  # Add subsequent layers ------------------------------------------------------
  for (j in 2:m) {
    p <- p +
      annotate(geom = "raster", x = plot.df[, 1], y = plot.df[, 2],
               alpha = alpha * plot.df[, class.ind[j]], fill = fill.col[j])
      # geom_line(data = decbound.df[[j]], aes(X1, X2))
  }

  # Add points and touch up remaining plot  ------------------------------------
  p <- p +
    geom_point(data = points.df, aes(X1, X2, col = Class)) +
    coord_cartesian(xlim = x, ylim = y) +
    guides(fill = FALSE, alpha = FALSE) +
    theme_bw()
  p
}

# iplot_decbound <- function(x, levels = NULL) {
#   classes <- as.factor(x$y)
#   if (!is.null(levels)) levels(classes) <- levels
#   else levels(classes) <- x$y.levels
#   xx <- c(min(x$X[, 1]), max(x$X[, 1]))
#   yy <- c(min(x$X[, 2]), max(x$X[, 2]))
#   decision.line <- boundarySolver(x$alpha, x$lambda, x$w, x$kernel, x$X,
#                                   xmin = xx[1] - 10, xmax = xx[2] + 10,
#                                   ymin = yy[1] - 100, ymax = yy[2] + 100)
#   plot.df <- data.frame(X = x$X, Observation = 1:nrow(x$X), Class = classes)
#   dec.df <- data.frame(x.dec = decision.line$x1,
#                        y.dec = decision.line$x2)
#
#   ggplot(data = plot.df, aes(x = X[, 1], y = X[, 2])) +
#     geom_point(aes(col = Class), size = 3) +
#     # geom_text(aes(label = Observation)) +
#     geom_line(data = dec.df, aes(x = x.dec, y = y.dec), col = "grey35", linetype = "longdash") +
#     coord_cartesian(xlim = xx, ylim = yy) +
#     labs(x = colnames(x$X)[1], y = colnames(x$X)[2]) +
#     theme_bw()
# }
#
# boundarySolver <- function(alpha, lambda, w, kernel, X, xmin, xmax,
#                            ymin, ymax) {
#   xlength <- 500
#   x1plot <- seq(xmin, xmax, length = xlength)
#   x2plot <- rep(NA, xlength)
#
#   objFn <- function(x2, x1) {
#     H.tilde <- ikernL(list(X), list(matrix(c(x1, x2), nrow = 1)),
#                       kernel = kernel)[[1]]
#     alpha + lambda * as.numeric(H.tilde %*% w)
#   }
#   i <- 1
#   while (i <= xlength) {
#     tmp <- uniroot(objFn, c(ymin, ymax), x1 = x1plot[i])
#     x2plot[i] <- tmp$root
#     if (x2plot[i] > ymax) break
#     i <- i + 1
#   }
#
#   list(x1 = x1plot, x2 = x2plot)
# }
