#' Plot intervals between consecutive samples
#'
#' \code{plotGaps()} draws the distribution of the intervals between
#' consecutive samples of each subject.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the intervals into groups. Each interval
#'   takes the value of its later sample. (Default: \code{NULL})
#'
#' @param layout \code{Character scalar}. Draws a histogram
#'   (\code{"histogram"}), density curves (\code{"density"}) or box plots
#'   (\code{"box"}). (Default: \code{"histogram"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_histogram]{geom_histogram}},
#'   \code{\link[ggplot2:geom_density]{geom_density}} or
#'   \code{\link[ggplot2:geom_boxplot]{geom_boxplot}}.
#'
#' @details
#' Intervals are in \code{time.col} units: days for \code{Date}, hours for
#' \code{POSIXct} times and visit steps for visits.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotGaps
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 2))]
#' plotGaps(tse, time.col = "time", group = "subject", bins = 30)
#'
#' # Intervals by sex
#' plotGaps(tse, time.col = "time", group = "subject", group.by = "sex",
#'     layout = "box")
NULL

#' @rdname plotGaps
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotGaps", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, group.by = NULL,
            layout = c("histogram", "density", "box"), ...) {
        layout <- match.arg(layout)
        .require_column(x, group)
        .check_longitudinal(x, time.col, group)
        .check_column(x, group.by)

        df <- .order_time(as.data.frame(colData(x), optional = TRUE),
            time.col, group)
        i <- .next_in_subject(df[[group]])
        if (length(i) == 0L) {
            stop("No subject has more than one sample.", call. = FALSE)
        }
        df <- data.frame(interval = diff(.time_numeric(df[[time.col]]))[i],
            df[i + 1L, group.by, drop = FALSE], check.names = FALSE)
        .plot_gaps(df, time.col, group.by, layout, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_histogram geom_density
#'   geom_boxplot labs
.plot_gaps <- function(df, time.col, group.by, layout, ...) {
    by <- !is.null(group.by)
    lab <- paste("interval in", time.col)
    p <- switch(layout,
        box = ggplot(df, aes(x = if (by) .data[[group.by]] else "",
            y = .data$interval)) + geom_boxplot(...) +
            labs(x = group.by, y = lab),
        histogram = ggplot(df, aes(x = .data$interval,
            fill = if (by) .data[[group.by]])) + geom_histogram(...) +
            labs(x = lab, fill = group.by),
        density = ggplot(df, aes(x = .data$interval,
            colour = if (by) .data[[group.by]])) + geom_density(...) +
            labs(x = lab, colour = group.by))
    p + .scale_by(df, group.by, c("colour", "fill")) + .theme_tempoPlot()
}
