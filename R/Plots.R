################################################################################
#
#   iprobit: Binary and Multinomial Probit Regression with I-priors
#   Copyright (C) 2017  Haziq Jamil
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

#' @export
plot.iprobitMod <- function(x, niter.plot = NULL, levels = NULL, ...) {
  iplot_fitted(x)
}

#' @export
iplot_fitted <- function(object) {
  list2env(object, environment())
  list2env(ipriorKernel, environment())
  list2env(model, environment())

  probs <- fitted(object)$prob
  if (isNystrom(ipriorKernel)) probs <- probs[order(Nystrom$Nys.samp), ]
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

#' @export
iplot_lb <- function(x, niter.plot = NULL, lab.pos = c("up", "down")) {
  if (x$niter < 2) stop("Nothing to plot.")

  lab.pos <- match.arg(lab.pos, c("up", "down"))
  if (lab.pos == "up") lab.pos <- -0.5
  else lab.pos <- 1.5

  lb.original <- x$lower.bound
  if (is.null(niter.plot)) niter.plot <- c(1, length(lb.original))
  else if (length(niter.plot) == 1) niter.plot <- c(1, niter.plot)
  niter.plot <- niter.plot[1]:niter.plot[2]
  lb <- lb.original[niter.plot]
  plot.df <- data.frame(Iteration = niter.plot, lb = lb)
  time.per.iter <- x$time$time / x$niter
  if (time.per.iter < 0.001) time.per.iter <- 0.001

  ggplot(plot.df, aes(x = Iteration, y = lb, label = max(lb))) +
    geom_line(col = "grey60") +
    geom_point() +
    geom_hline(yintercept = max(lb.original), linetype = 2, col = "red") +
    scale_x_continuous(
      sec.axis = sec_axis(~ . * time.per.iter, name = "Time (seconds)"),
      breaks = scales::pretty_breaks(n = min(5, ifelse(x$niter == 2, 1, x$niter)))
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

  if (isNystrom(object)) {
    Nys.m <- object$ipriorKernel$Nystrom$m
    Nys.lab <- object$ipriorKernel$Nystrom$Nys.samp[1:Nys.m]
    Nys.df <- points.df[seq_len(Nys.m), ]
    Nys.df <- cbind(Nys.df, Nys.lab)
    p <- p +
      geom_point(data = Nys.df, aes(X1, X2), size = 1.8) +
      geom_point(data = Nys.df, aes(X1, X2, col = Class), size = 1)
  }

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
    # annotate(geom = "raster", x = plot.df[, 1], y = plot.df[, 2],
    #          alpha = 0.6 * plot.df[, 3], fill = "#F8766D") +
    # annotate(geom = "raster", x = plot.df[, 1], y = plot.df[, 2],
    #          alpha = 0.6 * plot.df[, 4], fill = "#00BFC4") +
    geom_point(data = points.df, aes(X1, X2, col = Class)) +
               # col = "black", shape = 21, stroke = 0.8) +
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

#' @export
iplot_error <- function(x, niter.plot = NULL) {
  if (x$niter < 2) stop("Nothing to plot.")

  if (is.null(niter.plot)) niter.plot <- c(1, length(x$error))
  else if (length(niter.plot) == 1) niter.plot <- c(1, niter.plot)
  niter.plot <- niter.plot[1]:niter.plot[2]
  plot.df <- data.frame(Iteration = niter.plot,
                        error     = x$error[niter.plot] / 100,
                        brier     = x$brier[niter.plot])
  time.per.iter <- x$time$time / x$niter
  if (time.per.iter < 0.001) time.per.iter <- 0.001
  plot.df <- reshape2::melt(plot.df, id = "Iteration")
  last.value.df <- subset(plot.df, plot.df$Iteration == max(plot.df$Iteration))
  last.value.df$lab <- NA
  last.value.df$lab[1] <- decimal_place(last.value.df$value[1] * 100)
  last.value.df$lab[1] <- paste0(last.value.df$lab[1], "%")
  last.value.df$lab[2] <- decimal_place(last.value.df$value[2], 3)

  ggplot(plot.df, aes(x = Iteration, y = value, col = variable)) +
    geom_line(alpha = 0.5) +
    directlabels::geom_dl(aes(label = variable), method = ("lines2")) +
    geom_point(size = 0.9) +
    geom_text(data = last.value.df, aes(label = lab, y = value + 0.002),
              vjust = 0, hjust = 0.5) +
    scale_x_continuous(
      sec.axis = sec_axis(~ . * time.per.iter, name = "Time (seconds)"),
      breaks = scales::pretty_breaks(n = min(5, ifelse(x$niter == 2, 1, x$niter)))
    ) +
    scale_y_continuous(
      name = "Training error rate",
      labels = scales::percent,
      sec.axis = sec_axis(~ ., name = "Brier score")
    ) +
    theme_bw() +
    theme(legend.position = "none")
}

iplot_lb_and_error <- function(x, niter.plot = 10, lab.pos = c("up", "down")) {
  p1 <- iplot_lb(mod, niter.plot, lab.pos) +
    scale_y_continuous(sec.axis = dup_axis())
    # theme(axis.title.x=element_blank(),
    #       axis.text.x=element_blank(),
    #       axis.ticks.x=element_blank())
  suppressMessages(
    p2 <- iplot_error(mod, niter.plot) +
      scale_x_continuous(
        breaks = scales::pretty_breaks(n = min(5, ifelse(x$niter == 2, 1, x$niter)))
      )
  )

  cowplot::plot_grid(p1, p2, nrow = 2)
}

# tmp <- cowplot::plot_grid(iplot_lb(mod), NULL, rel_widths = c(1, 0.06))
# cowplot::plot_grid(tmp, iplot_error(mod), nrow = 2)

# iplot_prob <- function(x, covariate = 1, levels = NULL) {
#   if (!is.numeric(covariate)) {
#     covariate <- grep(covariate, colnames(mod$X))
#   }
#   classes <- as.factor(x$y)
#   if (!is.null(levels)) levels(classes) <- levels
#   else levels(classes) <- x$y.levels
#   mean.X <- apply(x$X, 2, mean)
#   X.new <- x$X
#   X.new[, -covariate] <- mean.X[-covariate]
#   p.hat <- predict(x, X.new)$prob
#   p.hat.lower <- predict(x, X.new, "lower")$prob
#   p.hat.upper <- predict(x, X.new, "upper")$prob
#   plot.df <- data.frame(x = x$X[, covariate], y = x$y, prob = p.hat,
#                         low = p.hat.lower, high = p.hat.upper,
#                         Class = classes)
#   x.axis.lab <- colnames(x$X)[covariate]
#
#   ggplot(plot.df, aes(x = x, y = y)) +
#     geom_line(aes(x = x, y = prob)) +
#     geom_point(aes(col = Class)) +
#     geom_ribbon(aes(x = x, ymin = low, ymax = high), fill = "grey70", alpha = 0.5) +
#     labs(x = x.axis.lab, y = "Probability") +
#     theme_bw()
# }
