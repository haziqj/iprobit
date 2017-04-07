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