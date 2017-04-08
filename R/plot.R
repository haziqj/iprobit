#' Plot for \code{ipriorProbit} objects.
#'
#' @param x
#' @param niter.plot
#' @param levels
#' @param ...
#'
#' @return
#' @export
plot.ipriorProbit <- function(x, niter.plot = NULL, levels = NULL, ...) {
  fitted.plot <- iplot_fitted(x, levels = levels)
}

iplot_fitted <- function(x, levels = NULL) {
  classes <- as.factor(x$y)
  if (!is.null(levels)) levels(classes) <- levels
  else levels(classes) <- x$y.levels
  plot.df <- data.frame(Observation = 1:length(x$ystar), p.hat = fitted(x)$prob,
                        Class = classes)

  ggplot(plot.df, aes(x = Observation, y = p.hat, col = Class)) +
    geom_point() +
    labs(y = "Fitted probabilities") +
    theme_bw()
}

iplot_lb <- function(x, niter.plot = NULL) {
  lb.original <- x$lower.bound[!is.na(x$lower.bound)]
  if (is.null(niter.plot)) niter.plot <- c(1, length(lb.original))
  else if (length(niter.plot) == 1) niter.plot <- c(1, niter.plot)
  niter.plot <- niter.plot[1]:niter.plot[2]
  lb <- lb.original[niter.plot]
  plot.df <- data.frame(Iteration = niter.plot, lb = lb)
  time.per.iter <- as.numeric(x$time) / (x$niter - 1)
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
             vjust = -0.5, label = round(max(lb.original), 2), size = 3.7) +
    labs(y = "Variational lower bound") +
    theme_bw()
}

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

iplot_decbound <- function(x, levels = NULL) {
  classes <- as.factor(x$y)
  if (!is.null(levels)) levels(classes) <- levels
  else levels(classes) <- x$y.levels
  plot.df <- data.frame(X = x$X, Class = classes)

  ggplot(data = plot.df, aes(x = X[, 1], y = X[, 2], col = Class)) +
    geom_point(size = 3) +
    labs(x = colnames(x$X)[1], y = colnames(x$X)[2]) +
    theme_bw()
}






