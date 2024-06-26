ggInteract_2D_test <-
function (gbm.object, dat.rank = NULL, n.plots = NULL, x = 1, 
          y = 2, col.gradient = c("lightskyblue", "#142c42"), palette = NULL, 
          show.axis = F, show.dot = F, contour = T, label.contour = F, 
          col.contour = "cadetblue2", col.dot = "grey10", alpha.dot = 0.9, 
          cex.dot = 0.5, pred.means = NULL, x.label = NULL, y.label = NULL, 
          legend = TRUE, z.label = "Fitted value", x.range = NULL, 
          nrow = NULL, ncol = NULL, y.range = NULL, z.range = NULL, 
          main = "", smooth = "none", ...) {
  gbm.call <- gbm.object$gbm.call
  pred.names <- gbm.call$predictor.names
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  gbm.y <- gbm.call$gbm.y
  family = gbm.call$family
  have.factor <- FALSE
  if (is.null(dat.rank)) {
    if (is.character(x)) {
      x <- match(x, pred.names)
    }
    if (is.character(y)) {
      y <- match(y, pred.names)
    }
    x.name <- gbm.call$predictor.names[x]
    if (is.null(x.label)) {
      x.label <- gbm.call$predictor.names[x]
    }
    y.name <- gbm.call$predictor.names[y]
    if (is.null(y.label)) {
      y.label <- gbm.call$predictor.names[y]
    }
    data <- gbm.call$dataframe[, gbm.x, drop = FALSE]
    n.trees <- gbm.call$best.trees
    if (is.vector(data[, x])) {
      if (is.null(x.range)) {
        x.var <- seq(min(data[, x], na.rm = T), max(data[, 
                                                         x], na.rm = T), length = 50)
      }
      else {
        x.var <- seq(x.range[1], x.range[2], length = 50)
      }
    }
    else {
      x.var <- names(table(data[, x]))
      have.factor <- TRUE
    }
    if (is.vector(data[, y])) {
      if (is.null(y.range)) {
        y.var <- seq(min(data[, y], na.rm = T), max(data[, 
                                                         y], na.rm = T), length = 50)
      }
      else {
        y.var <- seq(y.range[1], y.range[2], length = 50)
      }
    }
    else {
      y.var <- names(table(data[, y]))
      if (have.factor) {
        stop("at least one marginal predictor must be a vector!")
      }
      else {
        have.factor <- TRUE
      }
    }
    pred.frame <- expand.grid(list(x.var, y.var))
    names(pred.frame) <- c(x.name, y.name)
    pred.rows <- nrow(pred.frame)
    if (have.factor) {
      if (is.factor(pred.frame[, 2])) {
        pred.frame <- pred.frame[, c(2, 1)]
        x.var <- y.var
      }
    }
    j <- 3
    for (i in 1:n.preds) {
      if (i != x & i != y) {
        if (is.vector(data[, i])) {
          m <- match(pred.names[i], names(pred.means))
          if (is.na(m)) {
            pred.frame[, j] <- mean(data[, i], na.rm = T)
          }
          else pred.frame[, j] <- pred.means[m]
        }
        if (is.factor(data[, i])) {
          m <- match(pred.names[i], names(pred.means))
          temp.table <- table(data[, i])
          if (is.na(m)) {
            pred.frame[, j] <- rep(names(temp.table)[2], 
                                   pred.rows)
          }
          else {
            pred.frame[, j] <- pred.means[m]
          }
          pred.frame[, j] <- factor(pred.frame[, j], 
                                    levels = names(temp.table))
        }
        names(pred.frame)[j] <- pred.names[i]
        j <- j + 1
      }
    }
    prediction <- gbm::predict.gbm(gbm.object, pred.frame, 
                                   n.trees = n.trees, type = "response")
    if (smooth == "model") {
      pred.glm <- glm(prediction ~ ns(pred.frame[, 1], 
                                      df = 8) * ns(pred.frame[, 2], df = 8), data = pred.frame, 
                      family = poisson)
      prediction <- fitted(pred.glm)
    }
    max.pred <- max(prediction)
    cat("maximum value = ", round(max.pred, 2), "\n")
    if (is.null(z.range)) {
      if (family == "bernoulli") {
        z.range <- c(0, 1)
      }
      else if (family == "poisson") {
        z.range <- c(0, max.pred * 1.1)
      }
      else {
        z.min <- min(data[, y], na.rm = T)
        z.max <- max(data[, y], na.rm = T)
        z.delta <- z.max - z.min
        z.range <- c(z.min - (1.1 * z.delta), z.max + 
                       (1.1 * z.delta))
      }
    }
    if (have.factor == FALSE) {
      pred.matrix <- matrix(prediction, ncol = 50, nrow = 50)
      if (smooth == "average") {
        pred.matrix.smooth <- pred.matrix
        for (i in 2:49) {
          for (j in 2:49) {
            pred.matrix.smooth[i, j] <- mean(pred.matrix[c((i - 
                                                              1):(i + 1)), c((j - 1):(j + 1))])
          }
        }
        pred.matrix <- pred.matrix.smooth
      }
      rownames(pred.matrix) <- x.var
      colnames(pred.matrix) <- y.var
      ggdata <- melt(pred.matrix)
      pred.data1 <- gbm.call$dataframe[, gbm.call$gbm.x[x]]
      pred.data2 <- gbm.call$dataframe[, gbm.call$gbm.x[y]]
      df_point <- data.frame(pred1 = pred.data1, pred2 = pred.data2, 
                             value = gbm.object$fitted)
      p2D <- ggplot(data = ggdata, aes(x = Var1, y = Var2, 
                                       z = value)) + geom_tile(aes(fill = value), show.legend = legend) + 
        scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 
                                                                             0)) + labs(title = main, x = x.label, y = y.label)
      if (is.null(palette)) {
        p2D <- p2D + scale_fill_gradient(name = z.label, 
                                         limits = z.range, low = col.gradient[1], high = col.gradient[2])
      }
      else {
        p2D <- p2D + scale_fill_distiller(name = z.label, 
                                          limits = z.range, palette = palette, direction = 1)
      }
      if (show.dot == T) {
        p2D <- p2D + geom_point(data = df_point, aes(x = pred1, 
                                                     y = pred2), color = col.dot, alpha = alpha.dot, 
                                size = cex.dot)
      }
      if (show.axis == T) {
        p2D <- p2D + theme(axis.line.x = element_line(color = "gray28", 
                                                      size = 0.3), axis.line.y = element_line(color = "gray28", 
                                                                                              size = 0.3),
                           axis.title.x = element_text(size = 15), 
                           axis.title.y = element_text(size = 15),
                           axis.text.x = element_text(size = 12),
                           axis.text.y = element_text(size = 12))
      }
      if (contour == T) {
        p2D <- p2D + stat_contour(colour = col.contour,linewidth=1)
      }
      if ((contour == T) & (label.contour == T)) {
        p2D <- p2D + stat_contour(colour = col.contour) + 
          geom_dl(aes(label = ..level..), colour = col.contour, 
                  method = list("far.from.others.borders", 
                                "calc.boxes", "enlarge.box", hjust = 1, 
                                vjust = 1, box.color = "NA", fill = "transparent", 
                                "draw.rects"), stat = "contour")
      }
      p2D
    }
    else {
      ggdata <- data.frame(factor = pred.frame[, 1], x = pred.frame[, 
                                                                    2], y = prediction)
      p2DFactor <- ggplot(ggdata, aes(x = x, y = y, colour = factor)) + 
        geom_line(aes(linetype = factor)) + labs(title = main, 
                                                 x = y.label, y = z.label) + ylim(z.range) + theme_bw() + 
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
              legend.title = element_blank(),
              axis.title.x = element_text(size = 15), 
              axis.title.y = element_text(size = 15),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12))
      p2DFactor
    }
  }
  else {
    ggInteract2D.plots <- function(gbm.object) {
      ggInteract2D <- list()
      gbm.call <- gbm.object$gbm.call
      pred.names <- gbm.call$predictor.names
      gbm.x <- gbm.call$gbm.x
      n.preds <- length(gbm.x)
      gbm.y <- gbm.call$gbm.y
      family = gbm.call$family
      if (is.null(n.plots)) {
        n.plots <- nrow(dat.rank)
      }
      for (a in 1:n.plots) {
        have.factor <- FALSE
        x <- dat.rank[a, 1]
        y <- dat.rank[a, 3]
        x.name <- gbm.call$predictor.names[x]
        x.label <- gbm.call$predictor.names[x]
        y.name <- gbm.call$predictor.names[y]
        y.label <- gbm.call$predictor.names[y]
        data <- gbm.call$dataframe[, gbm.x, drop = FALSE]
        n.trees <- gbm.call$best.trees
        if (is.vector(data[, x])) {
          if (is.null(x.range)) {
            x.var <- seq(min(data[, x], na.rm = T), max(data[, 
                                                             x], na.rm = T), length = 50)
          }
          else {
            x.var <- seq(x.range[1], x.range[2], length = 50)
          }
        }
        else {
          x.var <- names(table(data[, x]))
          have.factor <- TRUE
        }
        if (is.vector(data[, y])) {
          if (is.null(y.range)) {
            y.var <- seq(min(data[, y], na.rm = T), max(data[, 
                                                             y], na.rm = T), length = 50)
          }
          else {
            y.var <- seq(y.range[1], y.range[2], length = 50)
          }
        }
        else {
          y.var <- names(table(data[, y]))
          if (have.factor) {
            stop("at least one marginal predictor must be a vector!")
          }
          else {
            have.factor <- TRUE
          }
        }
        pred.frame <- expand.grid(list(x.var, y.var))
        names(pred.frame) <- c(x.name, y.name)
        pred.rows <- nrow(pred.frame)
        if (have.factor) {
          if (is.factor(pred.frame[, 2])) {
            pred.frame <- pred.frame[, c(2, 1)]
            x.var <- y.var
          }
        }
        j <- 3
        for (i in 1:n.preds) {
          if (i != x & i != y) {
            if (is.vector(data[, i])) {
              m <- match(pred.names[i], names(pred.means))
              if (is.na(m)) {
                pred.frame[, j] <- mean(data[, i], na.rm = T)
              }
              else pred.frame[, j] <- pred.means[m]
            }
            if (is.factor(data[, i])) {
              m <- match(pred.names[i], names(pred.means))
              temp.table <- table(data[, i])
              if (is.na(m)) {
                pred.frame[, j] <- rep(names(temp.table)[2], 
                                       pred.rows)
              }
              else {
                pred.frame[, j] <- pred.means[m]
              }
              pred.frame[, j] <- factor(pred.frame[, 
                                                   j], levels = names(temp.table))
            }
            names(pred.frame)[j] <- pred.names[i]
            j <- j + 1
          }
        }
        prediction <- gbm::predict.gbm(gbm.object, pred.frame, 
                                       n.trees = n.trees, type = "response")
        if (smooth == "model") {
          pred.glm <- glm(prediction ~ ns(pred.frame[, 
                                                     1], df = 8) * ns(pred.frame[, 2], df = 8), 
                          data = pred.frame, family = poisson)
          prediction <- fitted(pred.glm)
        }
        max.pred <- max(prediction)
        cat("maximum value = ", round(max.pred, 2), "\n")
        if (is.null(z.range)) {
          if (family == "bernoulli") {
            z.range <- c(0, 1)
          }
          else if (family == "poisson") {
            z.range <- c(0, max.pred * 1.1)
          }
          else {
            z.min <- min(data[, y], na.rm = T)
            z.max <- max(data[, y], na.rm = T)
            z.delta <- z.max - z.min
            z.range <- c(z.min - (1.1 * z.delta), z.max + 
                           (1.1 * z.delta))
          }
        }
        if (have.factor == FALSE) {
          pred.matrix <- matrix(prediction, ncol = 50, 
                                nrow = 50)
          if (smooth == "average") {
            pred.matrix.smooth <- pred.matrix
            for (i in 2:49) {
              for (j in 2:49) {
                pred.matrix.smooth[i, j] <- mean(pred.matrix[c((i - 
                                                                  1):(i + 1)), c((j - 1):(j + 1))])
              }
            }
            pred.matrix <- pred.matrix.smooth
          }
          rownames(pred.matrix) <- x.var
          colnames(pred.matrix) <- y.var
          ggdata <- melt(pred.matrix)
          pred.data1 <- gbm.call$dataframe[, gbm.call$gbm.x[x]]
          pred.data2 <- gbm.call$dataframe[, gbm.call$gbm.x[y]]
          df_point <- data.frame(pred1 = pred.data1, 
                                 pred2 = pred.data2, value = gbm.object$fitted)
          ggInteract2D[[a]] <- ggplot(data = ggdata, 
                                      aes(x = Var1, y = Var2, z = value)) + geom_tile(aes(fill = value), 
                                                                                      show.legend = legend) + scale_y_continuous(expand = c(0, 
                                                                                                                                            0)) + scale_x_continuous(expand = c(0, 0)) + 
            labs(x = x.label, y = y.label)
          if (is.null(palette)) {
            ggInteract2D[[a]] <- ggInteract2D[[a]] + 
              scale_fill_gradient(name = z.label, limits = z.range, 
                                  low = col.gradient[1], high = col.gradient[2])
          }
          else {
            ggInteract2D[[a]] <- ggInteract2D[[a]] + 
              scale_fill_distiller(name = z.label, limits = z.range, 
                                   palette = palette, direction = 1)
          }
          if (show.dot == T) {
            ggInteract2D[[a]] <- ggInteract2D[[a]] + 
              geom_point(data = df_point, aes(x = pred1, 
                                              y = pred2), color = col.dot, alpha = alpha.dot, 
                         size = cex.dot)
          }
          if (show.axis == T) {
            ggInteract2D[[a]] <- ggInteract2D[[a]] + 
              theme(axis.line.x = element_line(color = "gray28", 
                                               size = 0.3), axis.line.y = element_line(color = "gray28", 
                                                                                       size = 0.3),
                    axis.title.x = element_text(size = 15), 
                    axis.title.y = element_text(size = 15),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12))
          }
          if (contour == T) {
            ggInteract2D[[a]] <- ggInteract2D[[a]] + 
              stat_contour(colour = col.contour, linewidth=1)
          }
          if ((contour == T) & (label.contour == T)) {
            if (n.plots > 20) {
              ggInteract2D[[a]] <- ggInteract2D[[a]] + 
                stat_contour(colour = col.contour)
            }
            else {
              ggInteract2D[[a]] <- ggInteract2D[[a]] + 
                stat_contour(colour = col.contour) + 
                geom_dl(aes(label = ..level..), colour = col.contour, 
                        method = list("far.from.others.borders", 
                                      "calc.boxes", "enlarge.box", hjust = 1, 
                                      vjust = 1, box.color = "NA", fill = "transparent", 
                                      "draw.rects"), stat = "contour")
            }
          }
        }
        else {
          ggdata <- data.frame(factor = pred.frame[, 
                                                   1], x = pred.frame[, 2], y = prediction)
          ggInteract2D[[a]] <- ggplot(ggdata, aes(x = x, 
                                                  y = y, colour = factor)) + geom_line(aes(linetype = factor)) + 
            labs(x = y.label, y = z.label) + ylim(z.range) + 
            theme_bw() + theme(panel.grid.minor = element_blank(), 
                               panel.grid.major = element_blank(), legend.title = element_blank(),
                               axis.title.x = element_text(size = 15), 
                               axis.title.y = element_text(size = 15),
                               axis.text.x = element_text(size = 12),
                               axis.text.y = element_text(size = 12))
        }
      }
      list(ggInteract2D = ggInteract2D)
    }
    plot <- ggInteract2D.plots(gbm.object)
    do.call(grid.arrange, c(plot$ggInteract2D, list(nrow = nrow, 
                                                    ncol = ncol)))
  }
}
