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
iplot_lb <- function(x, niter.plot = NULL, lab.pos = c("up", "down"), ...) {
  lae.check <- FALSE
  extra.opt <- list(...)
  if (isTRUE(extra.opt$lb.and.error)) {
    lae.check <- TRUE
    x.major <- extra.opt$x.major
    x.minor <- extra.opt$x.minor
  }

  if (x$niter < 2) stop("Nothing to plot.")

  lab.pos <- match.arg(lab.pos, c("up", "down"))
  if (lab.pos == "up") lab.pos <- -0.5
  else lab.pos <- 1.5

  lb.original <- x$lower.bound
  if (missing(niter.plot)) niter.plot <- seq_along(lb.original)
  else if (length(niter.plot) == 1) niter.plot <- c(1, niter.plot)
  lb <- lb.original[niter.plot]
  plot.df <- data.frame(Iteration = niter.plot, lb = lb)
  time.per.iter <- x$time$time / x$niter
  if (time.per.iter < 0.001) time.per.iter <- 0.001
  if (isTRUE(lae.check)) plot.df$Iteration <- plot.df$Iteration * time.per.iter

  p <- ggplot(plot.df, aes(x = Iteration, y = lb, label = max(lb))) +
    geom_line(col = "grey60") +
    geom_point() +
    geom_hline(yintercept = max(lb.original), linetype = 2, col = "red") +
    geom_text(aes(label = ifelse(lb == max(lb), round(max(lb), 2), "")),
              vjust = 1.5, size = 3.7) +
    annotate("text", col = "red", x = min(plot.df$Iteration), y = max(lb.original),
             vjust = lab.pos, label = round(max(lb.original), 2), size = 3.7) +
    labs(y = "Variational lower bound") +
    theme_bw()

  if (isTRUE(lae.check)) {
    p <- p +
      labs(x = "Time (seconds)") +
      scale_y_continuous(sec.axis = dup_axis()) +
      scale_x_continuous(
        position = "top",
        breaks = x.major * time.per.iter,
        minor_breaks = x.minor * time.per.iter
      ) +
      coord_cartesian(xlim = c(min(niter.plot) * time.per.iter,
                               (max(niter.plot) + 0.5) * time.per.iter))
  } else {
    p <- p +
      scale_x_continuous(
        sec.axis = sec_axis(~ . * time.per.iter, name = "Time (seconds)"),
        breaks = scales::pretty_breaks(n = min(5, ifelse(x$niter == 2, 1, x$niter)))
      )
  }

  p
}

prepare_point_range <- function(object, X.var = c(1, 2), grid.len = 10,
                                plot.test = TRUE) {
  # Helper function for iplot_predict() and iplot_dec_bound()
  X <- lapply(object$ipriorKernel$Xl, as.matrix)  # make all X a matrix
  col.info <- sapply(X, ncol)  # no. of cols of each matrix in the list
  p <- sum(col.info)  # How many X columns in total?
  Xl.min <- lapply(X, function(x) apply(x, 2, min))  # get min X values
  Xl.max <- lapply(X, function(x) apply(x, 2, max))  # get max X values
  Xl.minmax <- mapply(rbind, Xl.min, Xl.max, SIMPLIFY = FALSE)  # put together

  # Create grid for X values ---------------------------------------------------
  xx <- list(NULL)
  for (i in seq_along(Xl.minmax)) {
    tmp <- matrix(0, nrow = grid.len, ncol = ncol(Xl.minmax[[i]]))
    for (j in seq_len(ncol(Xl.minmax[[i]]))) {
      mm <- as.numeric(Xl.minmax[[i]][, j])
      tmp[, j] <- seq(from = mm[1] - 1, to = mm[2] + 1, length.out = grid.len)
    }
    xx[[i]] <- tmp
  }
  xx.all <- do.call(cbind, xx)  # combine into one big matrix
  xx.all.l <- split(xx.all, rep(1:ncol(xx.all), each = nrow(xx.all)))  # split into list
  xx.all <- expand.grid(xx.all.l)  # expand the grid

  col.ind2 <- cumsum(col.info)
  col.ind1 <- col.ind2 - col.info + 1
  for (i in seq_along(xx)) {
    xx[[i]] <- xx.all[, col.ind1[i]:col.ind2[i]]
  }

  # Build the plotting data frame (predicted probabilities) --------------------
  prob <- predict_iprobit(object, xx, NULL)$prob  # predicted probabilities
  plot.df <- xx.all[, X.var]
  plot.df <- data.frame(plot.df, prob)
  colnames(plot.df)[1:2] <- c("X1", "X2")
  colnames(plot.df)[-(1:p)] <- seq_along(object$ipriorKernel$y.levels)

  # Build the plotting data frame (observed data points) -----------------------
  classes <- factor(object$ipriorKernel$y)
  levels(classes) <- object$ipriorKernel$y.levels
  XX <- do.call(cbind, X)
  points.df <- data.frame(XX[, X.var], classes, "")
  colnames(points.df) <- c("X1", "X2", "Class", "prob")

  # Test data frame for plotting points ----------------------------------------
  test.fit <- object$test
  if (!is.null(test.fit) & isTRUE(plot.test)) {
    y.test <- factor(object$ipriorKernel$y.test)
    levels(y.test) <- object$ipriorKernel$y.levels
    x.test <- object$ipriorKernel$Xl.test
    x.test <- do.call(cbind, x.test)  # combine into one big matrix
    prob.ind <- test.fit$prob
    for (j in seq_len(ncol(prob.ind))) {
      prob.ind[, j] <- y.test == levels(y.test)[j]
    }
    probs.test <- unlist(test.fit$prob)[unlist(prob.ind)]
    test.df <- data.frame(x.test[, X.var], y.test, iprior::dec_plac(probs.test, 2))
    colnames(test.df) <- c("X1", "X2", "Class", "prob")
    points.df <- rbind(points.df, test.df)
  }

  # Others ---------------------------------------------------------------------
  mm <- t(do.call(cbind, Xl.minmax))  # max and min for each variable (rows)
  if (!is.null(object$ipriorKernel$formula)) {
    # Fitted using formula
    xname <- colnames(plot.df)[1:2] <-
      attr(object$ipriorKernel$terms, "term.labels")[X.var]
  } else {
    xname <- object$ipriorKernel$xname
  }

  list(plot.df = plot.df, points.df = points.df, mm = mm, xname = xname[X.var])
}

#' @export
iplot_dec_bound <- function(object, X.var = c(1, 2), col = "grey35", size = 0.8,
                            grid.len = 50, ...) {
  list2env(prepare_point_range(object, X.var, grid.len, FALSE),
           envir = environment())
  m <- get_m(object)

  plot.df <- reshape2::melt(plot.df, id.vars = c("X1", "X2"))
  # plot.df$value[plot.df$value > 1 / m] <- 1
  # plot.df$value[plot.df$value <= 1 / m] <- 0
  # head(plot.df)

  ggplot() +
    geom_point(data = points.df, aes(X1, X2, col = Class)) +
    geom_contour(data = plot.df, aes(X1, X2, z = value, group = variable,
                                     size = "Decision\nboundary"),
                 bins = 2, col = col, ...) +
    coord_cartesian(xlim = mm[1, ], ylim = mm[2, ]) +
    scale_colour_manual(values = c(iprior::gg_col_hue(m), "grey30")) +
    scale_size_manual(values = size, name = NULL) +
    guides(col = guide_legend(order = 1)) +
    # guides(col = guide_legend(override.aes = list(linetype = c(0, 0, 1),
    #                                               shape = c(19, 19, NA)))) +
    theme_bw()
}

#' @export
iplot_predict <- function(object, X.var = c(1, 2), grid.len = 50,
                          dec.bound = TRUE, plot.test = TRUE) {
  list2env(prepare_point_range(object, X.var, grid.len, plot.test),
           envir = environment())

  if (is.iprobitMod_bin(object))
    p <- iplot_predict_bin(plot.df, points.df, mm[1, ], mm[2, ], 2, dec.bound)
  if (is.iprobitMod_mult(object))
    p <- iplot_predict_mult(plot.df, points.df, mm[1, ], mm[2, ], get_m(object),
                            dec.bound)

  # if (isNystrom(object)) {
  #   Nys.m <- object$ipriorKernel$Nystrom$m
  #   Nys.lab <- object$ipriorKernel$Nystrom$Nys.samp[1:Nys.m]
  #   Nys.df <- points.df[seq_len(Nys.m), ]
  #   Nys.df <- cbind(Nys.df, Nys.lab)
  #   p <- p +
  #     geom_point(data = Nys.df, aes(X1, X2), size = 1.8) +
  #     geom_point(data = Nys.df, aes(X1, X2, col = Class), size = 1)
  # }

  yname <- ifelse(object$ipriorKernel$yname == "y", "Class",
                  object$ipriorKernel$yname)

  p + labs(x = xname[1], y = xname[2], col = yname)
}

iplot_predict_bin <- function(plot.df, points.df, x, y, m, dec.bound) {
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

iplot_predict_mult <- function(plot.df, points.df, x, y, m, dec.bound) {
  fill.col <- iprior::gg_col_hue(m)
  alpha <- 0.6
  class.ind <- sort(seq(from = ncol(plot.df), by = -1, length = m))

  # Add first layer ------------------------------------------------------------
  p <- ggplot() +
    geom_raster(data = plot.df, aes(X1, X2, fill = fill.col[1],
                                    alpha = alpha * plot.df[, class.ind[1]])) +
    scale_alpha_continuous(range = c(0, 0.5))

  # Add subsequent layers ------------------------------------------------------
  for (j in 2:m) {
    p <- p +
      annotate(geom = "raster", x = plot.df[, 1], y = plot.df[, 2],
               alpha = alpha * plot.df[, class.ind[j]], fill = fill.col[j])
  }

  # Add points and decision boundary, and touch up remaining plot --------------
  if (isTRUE(dec.bound)) {
    decbound.df <- reshape2::melt(plot.df, id.vars = c("X1", "X2"))
    p <- p +
      geom_contour(data = decbound.df, aes(X1, X2, z = value, group = variable,
                                           col = variable,
                                           size = "Decision\nboundary"),
                   # binwidth = 0.5 + 1e-12,
                   bins = 2,
                   linetype = "dashed") +
      scale_size_manual(values = 0.8, name = NULL) +
      scale_color_discrete(name = " Class") +
      coord_cartesian(xlim = x, ylim = y) +
      guides(fill = FALSE, alpha = FALSE,
             size = guide_legend(override.aes = list(linetype = 2, col = "grey35")),
             col = guide_legend(order = 1, override.aes = list(linetype = 0))) +
      theme_bw() +
      theme(legend.key.width = unit(2, "line"))
  } else {
    p <- p +
      coord_cartesian(xlim = x, ylim = y) +
      guides(fill = FALSE, alpha = FALSE) +
      theme_bw()
  }

  # Add points -----------------------------------------------------------------
  p <- p +
    ggrepel::geom_label_repel(data = points.df, segment.colour = "grey25",
                              box.padding = 0.9, show.legend = FALSE,
                              aes(X1, X2, col = Class, label = prob)) +
    geom_point(data = points.df, aes(X1, X2, col = Class)) +
    geom_point(data = subset(points.df, points.df$prob != ""), aes(X1, X2),
               shape = 1, col = "grey25")

  p
}

#' @export
iplot_error <- function(x, niter.plot, plot.test = TRUE) {
  if (x$niter < 2) stop("Nothing to plot.")
  if (missing(niter.plot)) niter.plot <- seq_along(x$train.error)
  else if (length(niter.plot) == 1) niter.plot <- c(1, niter.plot)

  # Prepare plotting data frame ------------------------------------------------
  plot.df <- data.frame(Iteration = niter.plot,
                        error     = x$train.error[niter.plot] / 100,
                        brier     = x$train.brier[niter.plot],
                        Type      = "train")
  if (!is.null(x$test) & isTRUE(plot.test)) {
    plot.df <- rbind(plot.df, data.frame(
      Iteration = niter.plot,
      error     = x$test.error[niter.plot] / 100,
      brier     = x$test.brier[niter.plot],
      Type      = "test"
    ))
  }
  time.per.iter <- x$time$time / x$niter
  if (time.per.iter < 0.001) time.per.iter <- 0.001
  plot.df <- reshape2::melt(plot.df, id = c("Iteration", "Type"))

  # Prepare data frame of the last points for labelling ------------------------
  last.value.df <- subset(plot.df, plot.df$Iteration == max(plot.df$Iteration))
  last.value.df$lab <- NA
  ind.error <- last.value.df$variable == "error"
  ind.brier <- last.value.df$variable == "brier"
  last.value.df$lab[ind.error] <-
    iprior::dec_plac(last.value.df$value[ind.error] * 100)
  last.value.df$lab[ind.error] <- paste0(
    "(", last.value.df$Type[ind.error], ")\n", last.value.df$lab[ind.error], "%"
  )
  last.value.df$lab[ind.brier] <- paste0(
    "(", last.value.df$Type[ind.brier], ")\n",
    iprior::dec_plac(last.value.df$value[ind.brier], 3)
  )

  # The plot -------------------------------------------------------------------
  ggplot(plot.df, aes(x = Iteration, y = value, col = variable,
                      linetype = Type)) +
    geom_line(alpha = 0.9) +
    directlabels::geom_dl(aes(label = variable), method = "extreme.grid") +
    geom_point(size = 0.9) +
    ggrepel::geom_text_repel(data = last.value.df, nudge_x = 0.5,
                             segment.colour = NA, size = 3.7,
                             aes(label = lab, y = value)) +
    scale_x_continuous(
      sec.axis = sec_axis(~ . * time.per.iter, name = "Time (seconds)"),
      breaks = scales::pretty_breaks(n = min(5, ifelse(x$niter == 2, 1, x$niter)))
    ) +
    scale_y_continuous(
      name = "Misclassification rate",
      labels = scales::percent,
      sec.axis = sec_axis(~ ., name = "Brier score")
    ) +
    coord_cartesian(xlim = c(min(niter.plot), max(niter.plot) + 0.5)) +
    theme_bw() +
    theme(legend.position = "none")
}

iplot_lb_and_error <- function(x, niter.plot, lab.pos) {
  suppressMessages(
    p2 <- iplot_error(x, niter.plot) +
      scale_x_continuous(
        breaks = scales::pretty_breaks(n = min(5, ifelse(x$niter == 2, 1, x$niter)))
      )
  )

  tmp <- ggplot_build(p2)$layout$panel_ranges[[1]]
  x.major <- tmp$x.major_source
  x.minor <- tmp$x.minor_source

  if (missing(lab.pos)) lab.pos <- "down"

  p1 <- iplot_lb(x, niter.plot, lab.pos, lb.and.error = TRUE,
                 x.major = x.major, x.minor = x.minor)

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
