#' @export
#' @keywords public
#'
#' @title Plot predictions for a FastRet model
#'
#' @description
#' Creates scatter plots of measured vs. predicted retention times (RT) for a
#' FastRet Model (FRM). Supports plotting cross-validation (CV) predictions and
#' fitted predictions on the training set, as well as their adjusted variants
#' when the model has been adjusted via [adjust_frm()]. Coloring highlights
#' points within 1 minute of the identity line and simple outliers.
#'
#' @param frm An object of class `frm` as returned by [train_frm()].
#' @param type Plot type. One of:
#' - "scatter.cv": CV predictions for the training set
#' - "scatter.cv.adj": CV predictions for the adjustment set (requires `frm$adj`)
#' - "scatter.train": Model predictions for the training set
#' - "scatter.train.adj": Adjusted model predictions for the adjustment set (requires `frm$adj`)
#' @param trafo Transformation applied for display. One of:
#' - "identity": no transformation
#' - "log2": apply log2 transform to axes (metrics are computed on raw values)
#'
#' @return NULL, called for its side effect of plotting.
#'
#' @examples
#' frm <- read_rp_lasso_model_rds()
#' plot_frm(frm, type = "scatter.cv")
plot_frm <- function(frm = train_frm(verbose = 1),
                     type = "scatter.cv",
                     trafo = "identity"
) {
    # Check args
    type <- match.arg(type, c("scatter.cv", "scatter.cv.adj", "scatter.train", "scatter.train.adj"))
    plot_adj <- grepl("adj", type)
    if (plot_adj && is.null(frm$adj)) {
        fmt <- paste(sep = "\n",
            "type is `%s`, but the model has not been adjusted yet.",
            "See `?adjust_frm` for information on how to adjust existing models."
        )
        stop(sprintf(fmt, type))
    }
    trafo <- match.arg(trafo, c("identity", "log2"))
    dotrafo <- switch(trafo, "identity" = identity, "log2" = log2)

    # Prepare data for plotting
    title <- switch(type,
        "scatter.cv" = "CV predictions for training data",
        "scatter.cv.adj" = "CV predictions for adjustment data",
        "scatter.train" = "Model predictions for training data",
        "scatter.train.adj" = "Adjusted model predictions for adjustment data"
    )
    df <- if (plot_adj) frm$adj$df else frm$df
    x <- if (plot_adj) df$RT_ADJ else df$RT
    y <- switch(type,
        "scatter.cv" = frm$cv$preds,
        "scatter.cv.adj" = frm$adj$cv$preds,
        "scatter.train" = predict(frm, df, adjust = FALSE),
        "scatter.train.adj" = predict(frm, df, adjust = TRUE)
    )
    is_within_1min <- (y > x - 1 & y < x + 1)
    is_min_outlier <- (y < (min(x) - mean(x)))
    is_max_outlier <- (y > (max(x) + mean(x)))
    is_outlier <- is_min_outlier | is_max_outlier
    col <- rep("blue", length(y))
    col[is_within_1min] <- "green"
    col[is_outlier] <- "red"
    y[is_min_outlier] <- min(x) - mean(x)
    y[is_max_outlier] <- max(x) + mean(x)
    R <- cor(x, y, on_zero_sd = NA)
    MSE <- mean((x - y)^2)
    MAE <- mean(abs(x - y))
    x <- dotrafo(x) # (1)
    y <- dotrafo(y)
    # (1) Important: do the transformation after calculating the performance
    # measures, as we only want to SHOW the transformed data, but not use it for
    # calculations

    # Plot data
    plot(
        x = x, y = y,
        pch = 21, # filled circle with border
        bg = col, # fill color
        col = "black", # border color
        lwd = 0.6, # border width
        xlab = if (trafo == "identity") "RT (Measured)" else "log2(RT) (Measured)",
        ylab = if (trafo == "identity") "RT (Predicted)" else "log2(RT) (Predicted)"
    )
    title(title)
    abline(a = 0, b = 1, col = "black")
    legend(
        x = "bottomright",
        legend = c(
            "Identity",
            sprintf("|y - x| < 1 min (n = %d)", sum(is_within_1min)),
            sprintf("|y - x| > 1 min (n = %d)", sum(!is_within_1min)),
            sprintf("y notin [min(x) - mu; max(x) + mu] (n = %d)", sum(is_outlier))
        ),
        col = c("black", "black", "black", "black"),
        pt.bg = c(NA, "green", "blue", "red"),
        pch = c(NA, 21, 21, 21),
        lty = c(1, 0, 0, 0),
        bg = grDevices::adjustcolor("white", alpha.f = 0.66)
    )
    legend(
        x = "topleft", col = "black", lty = 0, pch = NA,
        legend = c(sprintf("R = %.2f", R), sprintf("MSE = %.2f", MSE), sprintf("MAE = %.2f", MAE)),
        bg = grDevices::adjustcolor("white", alpha.f = 0.66)
    )
}

#' @noRd
#'
#' @description
#' Create boxplots for up to 4 models to compare their performance measures
#'
#' @param models
#' A list of objects of class `frm` as returned by [train_frm()].
#'
#' @param ptype
#' A string specifying the plot type. Options are: "base" (default) and "ggplot2".
#'
plot_boxplot <- function(model = train_frm(), ptype = "base") {
    .data <- rlang::.data # (1)
    # (1) To avoid warnings about NSE in ggplot2 calls
    # (https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html)
    ptype <- match.arg(ptype, c("base", "ggplot2"))
    stats <- as.data.frame(collect(model$cv$stats))
    p <- if (ptype == "base") {
        opar <- par(mar = c(2, 2, 2, 0) + 0.5)
        on.exit(par(opar), add = TRUE)
        boxplot(stats, main = "CV Performance across folds", xlab = "", ylab = "", cex.axis = 0.8)
    } else {
        stats_lf <- reshape(
            stats,
            varying = names(stats), v.names = "Value", timevar = "Measure",
            times = names(stats), direction = "long"
        )
        p <- ggplot2::ggplot(
            data = stats_lf,
            mapping = ggplot2::aes(x = .data$Measure, y = .data$Value)
        )
        p <- p + ggplot2::geom_boxplot()
        p <- p + ggplot2::ggtitle("CV Performance across folds")
        p <- p + ggplot2::theme_minimal()
        p <- p + ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )
    }
}
