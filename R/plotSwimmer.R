#' Plot sampling timelines of subjects
#'
#' \code{plotSwimmer()} draws one row per subject with time on the x axis
#' and one point per sample. A lane joins the samples of each subject.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param colour.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that colours the points. (Default: \code{NULL})
#'
#' @param order.by \code{Character scalar}. Orders subjects by their first
#'   sampling time (\code{"first"}) or by number of samples, most first
#'   (\code{"n"}). Ties keep the order of \code{group}. The first subject is
#'   at the top. (Default: \code{"first"})
#'
#' @param max.gap \code{Numeric scalar}. Breaks a lane where consecutive
#'   samples of a subject are further apart than this, in \code{time.col}
#'   units. \code{NULL} never breaks lanes. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Numeric and \code{Date} times are drawn on a continuous axis. Visits are
#' drawn in level order, including levels without samples.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotSwimmer
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
#' plotSwimmer(tse, time.col = "time", group = "subject",
#'     colour.by = "bmi_group")
NULL

#' @rdname plotSwimmer
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotSwimmer", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, colour.by = NULL,
            order.by = c("first", "n"), max.gap = NULL, ...) {
        order.by <- match.arg(order.by)
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        # Errors on character and other unsupported time types.
        .time_type(x, time.col)
        .check_column(x, colour.by)

        cols <- unique(c(time.col, group, colour.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df <- .order_time(df, time.col, group)
        df <- .break_gaps(df, time.col, group, max.gap)

        # Rows are grouped by subject, so each run is one subject.
        runs <- rle(as.character(df[[group]]))
        first <- cumsum(runs$lengths) - runs$lengths + 1L
        key <- switch(order.by,
            first = .time_numeric(df[[time.col]])[first],
            n = -runs$lengths)
        df[[group]] <- factor(df[[group]],
            levels = rev(runs$values[order(key)]))

        .plot_swimmer(df, time.col, group, colour.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs
#'   scale_x_discrete
.plot_swimmer <- function(df, time.col, group, colour.by = NULL, ...) {
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data[[group]])) +
        geom_line(aes(group = interaction(.data[[group]], .data$.segment))) +
        geom_point(aes(colour = if (!is.null(colour.by)) .data[[colour.by]]),
            ...) +
        labs(x = time.col, y = group, colour = colour.by)
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p
}
