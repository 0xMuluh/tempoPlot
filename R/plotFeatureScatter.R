#' Plot two statistics of each feature with a fit
#'
#' \code{plotFeatureScatter()} draws one point per feature from two
#' statistics, computed from an assay or read from \code{rowData(x)}, and
#' optionally a fitted curve: Taylor's power law of mean and variance, or
#' the neutral model of occupancy and mean relative abundance.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param x.stat,y.stat \code{Character scalar}. Statistics on the axes:
#'   \code{"mean"}, \code{"variance"}, \code{"cv"} (coefficient of
#'   variation) or \code{"occupancy"} (share of samples with a value above
#'   zero) of \code{assay.type}, or a numeric column of \code{rowData(x)}.
#'
#' @param assay.type \code{Character scalar}. Specifies the assay the
#'   statistics are computed from. (Default: \code{"counts"})
#'
#' @param fit \code{Character scalar}. Draws no fit (\code{"none"}), a
#'   straight line on log-log axes (\code{"power"}, Taylor's law), or the
#'   neutral model (\code{"neutral"}), which needs \code{x.stat = "mean"}
#'   and \code{y.stat = "occupancy"} of counts. (Default: \code{"none"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the samples; statistics and fits are
#'   computed within each group. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' With \code{fit = "power"}, both axes are log10 and the legend gives the
#' slope. With \code{fit = "neutral"}, the mean is of relative abundances
#' and the curve is \eqn{1 - F(1/N; Nmp, Nm(1-p))}, with \eqn{F} the beta
#' distribution function, \eqn{N} the mean sample total, \eqn{p} the mean
#' relative abundance and the migration rate \eqn{m} fitted by least
#' squares. With a fit, features with a zero statistic on a log axis are
#' not drawn and are left out of the fit.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Taylor LR (1961). Aggregation, variance and the mean. \emph{Nature}
#' 189:732-735.
#'
#' Sloan WT, et al. (2006). Quantifying the roles of immigration and chance
#' in shaping prokaryote community structure. \emph{Environmental
#' Microbiology} 8(4):732-740.
#'
#' @name plotFeatureScatter
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006[, hitchip1006$time == 0]
#'
#' # Taylor's law of the genera
#' plotFeatureScatter(tse, x.stat = "mean", y.stat = "variance",
#'     fit = "power")
#'
#' # Occupancy against mean relative abundance in four artificial guts
#' data("SilvermanAGutData", package = "miaTime")
#' plotFeatureScatter(SilvermanAGutData, x.stat = "mean",
#'     y.stat = "occupancy", fit = "neutral", group.by = "Vessel")
NULL

#' @rdname plotFeatureScatter
#' @importFrom SummarizedExperiment colData assay assayNames rowData
#' @importFrom stats sd var
#' @importFrom methods setMethod
#' @export
setMethod("plotFeatureScatter", signature = c(x = "SummarizedExperiment"),
    function(x, x.stat, y.stat, assay.type = "counts",
            fit = c("none", "power", "neutral"), group.by = NULL, ...) {
        fit <- match.arg(fit)
        .check_column(x, group.by)
        if (!.is_one_name(assay.type) || !assay.type %in% assayNames(x)) {
            stop("'assay.type' must be one of assayNames(x).",
                call. = FALSE)
        }
        computed <- c("mean", "variance", "cv", "occupancy")
        for (arg in c("x.stat", "y.stat")) {
            v <- get(arg)
            if (!.is_one_name(v) || !(v %in% computed ||
                    (v %in% names(rowData(x)) &&
                    is.numeric(rowData(x)[[v]])))) {
                stop("'", arg, "' must be one of ",
                    paste(computed, collapse = ", "),
                    " or a numeric column of rowData(x).", call. = FALSE)
            }
        }
        if (fit == "neutral" && !(x.stat == "mean" &&
                y.stat == "occupancy")) {
            stop("fit = \"neutral\" needs x.stat = \"mean\" and ",
                "y.stat = \"occupancy\".", call. = FALSE)
        }
        mat <- as.matrix(assay(x, assay.type))
        by <- if (is.null(group.by)) rep("", ncol(x)) else
            as.character(colData(x)[[group.by]])
        res <- lapply(split(seq_len(ncol(x)), by), function(i) {
            m <- mat[, i, drop = FALSE]
            stat <- function(s) {
                if (!s %in% computed) return(rowData(x)[[s]])
                if (s == "mean" && fit == "neutral") {
                    return(rowMeans(prop.table(m, 2L), na.rm = TRUE))
                }
                switch(s, mean = rowMeans(m, na.rm = TRUE),
                    variance = apply(m, 1L, var, na.rm = TRUE),
                    cv = apply(m, 1L, sd, na.rm = TRUE) /
                        rowMeans(m, na.rm = TRUE),
                    occupancy = rowMeans(m > 0, na.rm = TRUE))
            }
            df <- data.frame(
                .feature = rownames(x) %||% as.character(seq_len(nrow(x))),
                .x = stat(x.stat),
                .y = stat(y.stat))
            curve <- switch(fit, none = NULL,
                power = .fit_power(df),
                neutral = .fit_neutral(df, mean(colSums(m, na.rm = TRUE))))
            list(points = df, curve = curve)
        })
        points <- .bind_groups(res, "points")
        curves <- .bind_groups(res, "curve")
        # Log axes cannot show zeros.
        if (fit != "none") {
            keep <- points$.x > 0 & (fit == "neutral" | points$.y > 0)
            points <- points[!is.na(keep) & keep, , drop = FALSE]
        }
        .plot_feature_scatter(points, curves, x.stat, y.stat, fit,
            group.by, ...)
    }
)

# Row-bind one element of each group's results, adding the group as .by.
.bind_groups <- function(res, what) {
    parts <- lapply(seq_along(res), function(i) {
        d <- res[[i]][[what]]
        if (is.null(d) || nrow(d) == 0L) return(NULL)
        cbind(.by = names(res)[i], d)
    })
    do.call(rbind, parts)
}

# Straight line on log10 axes; .label gives the slope.
#' @importFrom stats lm coef
.fit_power <- function(df) {
    ok <- is.finite(log10(df$.x)) & is.finite(log10(df$.y))
    if (sum(ok) < 3L) return(NULL)
    b <- coef(lm(log10(df$.y[ok]) ~ log10(df$.x[ok])))
    x <- 10^seq(min(log10(df$.x[ok])), max(log10(df$.x[ok])),
        length.out = 100L)
    data.frame(.x = x, .y = 10^(b[[1L]] + b[[2L]] * log10(x)),
        .label = sprintf("slope %.2f", b[[2L]]))
}

# Neutral model of Sloan et al. (2006) with the migration rate m fitted
# by least squares on log10(m).
#' @importFrom stats pbeta optimize
.fit_neutral <- function(df, n) {
    ok <- df$.x > 0 & df$.x < 1 & !is.na(df$.y)
    if (sum(ok) < 3L || !is.finite(n) || n <= 0) return(NULL)
    p <- df$.x[ok]
    occ <- df$.y[ok]
    model <- function(m, p) 1 - pbeta(1 / n, n * m * p, n * m * (1 - p))
    loss <- function(lm) sum((occ - model(10^lm, p))^2)
    lm <- optimize(loss, c(-8, 4))$minimum
    if (lm < -7.9 || lm > 3.9) {
        warning("The neutral model fit reached the bound of the migration ",
            "rate; the curve does not describe the data.", call. = FALSE)
    }
    m <- 10^lm
    x <- 10^seq(log10(min(p)), log10(max(p)), length.out = 100L)
    data.frame(.x = x, .y = model(m, x), .label = sprintf("m = %.2g", m))
}

#' @importFrom ggplot2 ggplot aes .data geom_point geom_line labs
#'   scale_x_log10 scale_y_log10
.plot_feature_scatter <- function(points, curves, x.stat, y.stat, fit,
        group.by, ...) {
    p <- ggplot(points, aes(x = .data$.x, y = .data$.y,
            colour = if (!is.null(group.by)) .data$.by)) +
        geom_point(...)
    if (!is.null(curves)) {
        curves$.curve <- if (is.null(group.by)) curves$.label else
            paste0(curves$.by, ": ", curves$.label)
        p <- p + geom_line(aes(x = .data$.x, y = .data$.y,
                group = .data$.by, linetype = .data$.curve),
            data = curves, inherit.aes = FALSE, colour = "grey20") +
            labs(linetype = "fit")
    }
    if (fit == "power") {
        p <- p + scale_y_log10()
    }
    if (fit %in% c("power", "neutral")) {
        p <- p + scale_x_log10()
    }
    xlab <- if (fit == "neutral") "mean relative abundance" else x.stat
    p + labs(x = xlab, y = y.stat, colour = group.by)
}
