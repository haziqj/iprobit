# myfun <- function() {
#   this.env <- environment()
#   environment(.lambdaExpand) <- this.env
#   list2env(mod, this.env)
#   list2env(BlockBstuff, this.env)
#   environment(BlockB)  <- this.env
#   list2env(model, this.env)
#   lambda <- 1:4
#   .lambdaExpand(env = this.env)
#   BlockB(1)
#   print(Sl)
# }

# gg_colour_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
#
# dat <- gen_circle(500, sd = 0.2)
# mod <- iprobit(dat$y, dat$X, kernel = "FBM", interactions = FALSE, maxit = 100)
#
# dat.test <- gen_circle(100, sd = 0.2)
# y.predict <- predict(mod, newdata = dat.test$X)$y
# mean(y.predict != dat.test$y) * 100
#
# p <- plot(dat)
# tmp <- cbind(apply(dat$X, 2, min), apply(dat$X, 2, max))
# x <- tmp[1, ]
# y <- tmp[2, ]
# xx <- seq(from = x[1] - 1, to = x[2] + 1, length.out = 100)
# yy <- seq(from = y[1] - 1, to = y[2] + 1, length.out = 100)
# plot.df <- expand.grid(xx, yy)
# prob <- predict(mod, newdata = plot.df)$prob
# plot.df <- cbind(plot.df, prob)
# colnames(plot.df) <- c("X1", "X2", "class1", "class2")
#
#
# ggplot() +
#   geom_raster(data = plot.df, aes(X1, X2, fill = class2 ), alpha = 0.5) +
#   scale_fill_gradient(low = "#F8766D", high = "#00BFC4", limits = c(0, 1)) +
#   geom_point(data = data.frame(X = dat$X, Class = dat$y), aes(X.1, X.2, col = Class)) +
#   coord_cartesian(xlim = x, ylim = y) +
#   theme_bw()
#
#
#
#
#
#
# dat <- gen_mixture(500, m = 4)
# mod <- iprobit_mult(dat$y, dat$X, kernel = "FBM", maxit = 100, common.RKHS.scale = FALSE)
#
# dat.test <- gen_mixture(100, m = 4)
# predict(mod, dat.test$X, dat.test$y)
#
# p <- plot(dat)
# tmp <- cbind(apply(dat$X, 2, min), apply(dat$X, 2, max))
# x <- tmp[1, ]
# y <- tmp[2, ]
# xx <- seq(from = x[1] - 1, to = x[2] + 1, length.out = 100)
# yy <- seq(from = y[1] - 1, to = y[2] + 1, length.out = 100)
# plot.df <- expand.grid(xx, yy)
# mod.pred <- predict(mod, plot.df)
# plot.df <- cbind(plot.df, mod.pred$prob, mod.pred$y)
# colnames(plot.df) <- c("X1", "X2", "class1", "class2", "class3", "class4", "y")
# plot.df$kol <- apply(plot.df[1:10, c(3, 4, 5, 6, 7)], 1,
#                      function(x) as.numeric(x[as.numeric(x[5])]))
#
# ggplot() +
#   geom_raster(data = plot.df, aes(X1, X2, fill = as.numeric(y), alpha = 0.5)) +
#   scale_fill_gradientn(colours = c("#FFFFFF", gg_colour_hue(4)), limits = c(0, 4)) +
#   guides(fill = FALSE, alpha = FALSE) +
#   geom_point(data = data.frame(X = dat$X, Class = dat$y), aes(X.1, X.2, col = Class)) +
#   coord_cartesian(xlim = x, ylim = y) +
#   theme_bw()
