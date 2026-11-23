#' Plot time intervals of features
#'
#' \code{plotIntervals()} draws time intervals of each feature as
#' horizontal bars, such as the intervals in which a feature differs
#' between groups.
#'
#' @param x a \code{data.frame} of intervals, one row per interval with
#'   columns \code{start} and \code{end}, or a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with
#'   such a table in its \code{metadata}.
#'
#' @param name \code{Character scalar}. Specifies the entry of
#'   \code{metadata(x)} holding the intervals.
#'
#' @param feature.col \code{Character scalar}. Specifies the column of
#'   feature names. (Default: \code{"feature"})
#'
#' @param colour.by \code{Character scalar}. Specifies the column that
#'   colours the intervals, such as the group with the higher values.
#'   (Default: \code{NULL})
#'
#' @param facet.by \code{Character scalar}. Specifies the column that
#'   splits the intervals into panels, such as the method that found them.
#'   (Default: \code{NULL})
#'
#' @param features \code{Character vector}. Features to plot. If
#'   \code{NULL}, plots the \code{top} features with the longest total
#'   time in intervals. (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features to plot when
#'   \code{features} is \code{NULL}. (Default: \code{20})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_segment]{geom_segment}}.
#'
#' @details
#' \code{start} and \code{end} are numeric, \code{Date}, \code{POSIXct} or
#' factors of visits with the same levels. An interval with \code{end}
#' equal to \code{start}, or missing, covers one time point and is drawn as
#' a point. Features are ordered by the start of their first interval, the
#' earliest at the top.
#'
#' The columns \code{feature}, \code{start}, \code{end} and
#' \code{dominant} match the significant intervals returned by MetaLonDA
#' (Metwally et al. 2018, doi:10.1186/s40168-018-0402-y); intervals from
#' SplinectomeR or from comparisons of fitted curves can be given in the
#' same form.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotIntervals
#'
#' @examples
#' # Illustrative significant intervals in the form returned by MetaLonDA:
#' # days in which each genus differs between two groups, and the group
#' # with the higher abundance
#' intervals <- data.frame(
#'     feature = c("Bacteroides", "Bacteroides", "Prevotella",
#'         "Akkermansia", "Bifidobacterium"),
#'     start = c(10, 60, 30, 0, 45),
#'     end = c(25, 90, 75, 20, 45),
#'     dominant = c("case", "control", "case", "control", "case"))
#' plotIntervals(intervals, colour.by = "dominant")
#'
#' # The same intervals stored in metadata(tse)
#' data("minimalgut", package = "miaTime")
#' tse <- minimalgut
#' S4Vectors::metadata(tse)$intervals <- intervals
#' plotIntervals(tse, name = "intervals", colour.by = "dominant")
NULL

#' @rdname plotIntervals
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotIntervals", signature = c(x = "SummarizedExperiment"),
    function(x, name, ...) {
        if (!.is_one_name(name) || !is.data.frame(metadata(x)[[name]])) {
            stop("'name' must name a table in metadata(x).", call. = FALSE)
        }
        plotIntervals(metadata(x)[[name]], ...)
    }
)

#' @rdname plotIntervals
#' @importFrom methods setMethod
#' @export
setMethod("plotIntervals", signature = c(x = "data.frame"),
    function(x, feature.col = "feature", colour.by = NULL,
            facet.by = NULL, features = NULL, top = 20, ...) {
        df <- .check_intervals(as.data.frame(x), feature.col, colour.by,
            facet.by)
        # Features with the longest total time in intervals.
        total <- tapply(.time_numeric(df$.end) - .time_numeric(df$.start),
            df$.feature, sum)
        mat <- matrix(total, dimnames = list(names(total), NULL))
        df <- df[df$.feature %in% .select_features(mat, features, top), ,
            drop = FALSE]
        first <- tapply(.time_numeric(df$.start), df$.feature, min)
        df$.feature <- factor(df$.feature,
            levels = rev(names(first)[order(first, names(first))]))
        .plot_intervals(df, colour.by, facet.by,
            if (is.factor(x$start)) levels(x$start), ...)
    }
)

# Check an interval table x and return it with columns .feature, .start
# and .end, in the class of start (visits as level codes). A missing end
# is the start.
.check_intervals <- function(x, feature.col, colour.by, facet.by) {
    .check_table_columns(x, feature.col = feature.col, null = FALSE)
    .check_table_columns(x, colour.by = colour.by, facet.by = facet.by)
    if (!all(c("start", "end") %in% names(x))) {
        stop("'x' must have columns 'start' and 'end'.", call. = FALSE)
    }
    start <- x$start
    end <- x$end
    type <- .time_type(start)
    same <- if (type == "visit") is.factor(end) &&
        identical(levels(end), levels(start)) else
        identical(.time_type(end), type)
    if (!same) {
        stop("'end' must have the type of 'start'",
            if (type == "visit") ", with the same levels", ".",
            call. = FALSE)
    }
    if (anyNA(start)) {
        stop("'start' must not contain NA values.", call. = FALSE)
    }
    end[is.na(end)] <- start[is.na(end)]
    x$.feature <- as.character(x[[feature.col]])
    x$.start <- .time_numeric(start)
    x$.end <- .time_numeric(end)
    if (any(x$.end < x$.start)) {
        stop("'end' must not be before 'start'.", call. = FALSE)
    }
    x$.start <- .time_restore(x$.start, start)
    x$.end <- .time_restore(x$.end, start)
    x
}

#' @importFrom ggplot2 ggplot aes .data geom_segment geom_point labs
#'   facet_wrap vars scale_x_continuous scale_y_discrete
#' @importFrom utils modifyList
.plot_intervals <- function(df, colour.by, facet.by, levels, ...) {
    # Without colour.by, all intervals take the first colour.
    fixed <- if (is.null(colour.by)) list(colour = .okabe_ito[[1L]])
    point <- df$.end == df$.start
    p <- ggplot(df, aes(y = .data$.feature,
            colour = if (!is.null(colour.by)) .data[[colour.by]])) +
        do.call(geom_segment, c(list(aes(x = .data$.start,
            xend = .data$.end, yend = .data$.feature),
            data = df[!point, , drop = FALSE]),
            modifyList(c(list(linewidth = 3), fixed), list(...)))) +
        do.call(geom_point, c(list(aes(x = .data$.start),
            data = df[point, , drop = FALSE], size = 3), fixed)) +
        # Bars and points are separate layers: keep the order of all
        # features.
        scale_y_discrete(drop = FALSE) +
        labs(x = "time", y = NULL, colour = colour.by) +
        .theme_tempoPlot()
    if (!is.null(colour.by)) {
        p <- p + .scale_by(df, colour.by, "colour")
    }
    if (!is.null(levels)) {
        p <- p + scale_x_continuous(breaks = seq_along(levels),
            labels = levels)
    }
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    p
}
