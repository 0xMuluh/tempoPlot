#' Plot longitudinal trends
#'
#' \code{plotTrend()} plots the trajectory of a feature, colData column or
#' ordination component over time for each subject, with optional group
#' summaries, boxplot paired layouts, and baseline subtraction.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param y \code{Character scalar}. Specifies the feature, column or
#'   ordination component to plot: a column of \code{colData(x)}, a row of
#'   \code{assay(x, assay.type)}, or \code{"<dimred>:<component>"}.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar} or \code{numeric scalar}.
#'   Specifies the assay to retrieve \code{y} from when \code{y} is a
#'   feature name. (Default: \code{NULL})
#'
#' @param layout \code{Character scalar}. Specifies the plot layout:
#'   \code{"line"} draws one line per subject; \code{"box"} draws boxes of
#'   \code{y} per time point. (Default: \code{"line"})
#'
#' @param show.individuals \code{Logical scalar}. When \code{TRUE}, draws
#'   subject lines and points. When \code{FALSE}, hides subject lines and
#'   points, leaving only boxes or the smooth curve. (Default: \code{TRUE})
#'
#' @param smooth \code{Character scalar}. Specifies a group summary on top:
#'   \code{"none"}, \code{"mean"} (mean at each time point), or
#'   \code{"loess"} (loess curve). (Default: \code{"none"})
#'
#' @param baseline \code{Logical scalar}. When \code{TRUE}, subtracts each
#'   subject's first value from their trajectory. (Default: \code{FALSE})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} to colour subjects and smooth curves.
#'   (Default: \code{NULL})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} to facet by. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Breaks a trajectory where consecutive
#'   samples of a subject are further apart than this, in \code{time.col}
#'   units. \code{NULL} never breaks trajectories. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to the main layer
#'   (\code{\link[ggplot2:geom_line]{geom_line}} for \code{layout = "line"}
#'   or \code{\link[ggplot2:geom_boxplot]{geom_boxplot}} for
#'   \code{layout = "box"}).
#'
#' @details
#' Numeric and \code{Date} times are drawn on a continuous axis. Visits are
#' drawn in level order, including levels without samples.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotTrend
#'
#' @examples
#' data("minimalgut", package = "miaTime")
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier")
#'
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", smooth = "loess")
#'
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", layout = "box")
NULL

#' @rdname plotTrend
#' @importFrom SummarizedExperiment colData
#' @importFrom stats ave
#' @importFrom methods setMethod
#' @export
setMethod("plotTrend", signature = c(x = "SummarizedExperiment"),
    function(x, y, time.col, group = NULL, assay.type = NULL,
            layout = c("line", "box"), show.individuals = TRUE,
            smooth = c("none", "mean", "loess"), baseline = FALSE,
            group.by = NULL, facet.by = NULL, max.gap = NULL, ...) {
        layout <- match.arg(layout)
        smooth <- match.arg(smooth)
        if (!(isTRUE(show.individuals) || isFALSE(show.individuals))) {
            stop("'show.individuals' must be TRUE or FALSE.", call. = FALSE)
        }
        if (!(isTRUE(baseline) || isFALSE(baseline))) {
            stop("'baseline' must be TRUE or FALSE.", call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        # Errors on character and other unsupported time types.
        .time_type(x, time.col)
        .check_column(x, group.by)
        .check_column(x, facet.by)

        cols <- unique(c(time.col, group, group.by, facet.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df$.y <- unname(.resolve_y(x, y, assay.type))
        df <- .order_time(df, time.col, group)
        df <- .break_gaps(df, time.col, group, max.gap)

        ids <- if (is.null(group)) rep(1L, nrow(df)) else df[[group]]
        if (baseline) {
            df$.y <- df$.y - ave(df$.y, ids, FUN = function(v) v[1L])
        }
        df$.line <- interaction(ids, df$.segment, drop = TRUE)

        .plot_trend(df, time.col, y, layout, show.individuals, smooth,
            group.by, facet.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point geom_boxplot
#'   geom_smooth stat_summary facet_wrap vars labs scale_x_discrete
.plot_trend <- function(df, time.col, y, layout, show.individuals, smooth,
        group.by = NULL, facet.by = NULL, ...) {
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data$.y,
        colour = if (!is.null(group.by)) .data[[group.by]]))
    if (layout == "box") {
        p <- p + geom_boxplot(aes(group = .data[[time.col]], colour = NULL),
            ...)
    }
    if (show.individuals) {
        lines <- aes(group = .data$.line)
        p <- p + if (layout == "line") geom_line(lines, ...) else
            geom_line(lines)
        p <- p + geom_point()
    }
    summary_group <- aes(group = if (is.null(group.by)) 1L else
        .data[[group.by]])
    if (smooth == "mean") {
        p <- p + stat_summary(summary_group, fun = mean, geom = "line")
    } else if (smooth == "loess") {
        p <- p + geom_smooth(summary_group, method = "loess",
            formula = y ~ x, se = FALSE)
    }
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p + labs(x = time.col, y = y, colour = group.by)
}
