#' Plot intervals between consecutive samples
#'
#' \code{plotGaps()} draws the distribution of the intervals between
#' consecutive samples of each subject.
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
#' Intervals are in \code{time.col} units: days for \code{Date} times and
#' visit steps for visits.
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
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        .check_column(x, group.by)

        cols <- unique(c(time.col, group, group.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df <- .order_time(df, time.col, group)
        n <- nrow(df)
        same <- df[[group]][-1L] == df[[group]][-n]
        df <- data.frame(
            interval = diff(.time_numeric(df[[time.col]]))[same],
            df[-1L, group.by, drop = FALSE][same, , drop = FALSE],
            check.names = FALSE)
        if (nrow(df) == 0L) {
            stop("No subject has more than one sample.", call. = FALSE)
        }
        .plot_gaps(df, time.col, group.by, layout, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_histogram geom_density
#'   geom_boxplot labs
.plot_gaps <- function(df, time.col, group.by, layout, ...) {
    lab <- paste("interval in", time.col)
    if (layout == "box") {
        return(ggplot(df, aes(x = if (is.null(group.by)) "" else
                .data[[group.by]], y = .data$interval)) +
            geom_boxplot(...) +
            labs(x = group.by, y = lab))
    }
    p <- ggplot(df, aes(x = .data$interval)) + labs(x = lab)
    if (layout == "histogram") {
        return(p + geom_histogram(aes(fill = if (!is.null(group.by))
            .data[[group.by]]), ...) + labs(fill = group.by))
    }
    p + geom_density(aes(colour = if (!is.null(group.by))
        .data[[group.by]]), ...) + labs(colour = group.by)
}
