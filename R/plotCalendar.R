#' Plot sample counts on a calendar
#'
#' \code{plotCalendar()} counts the samples taken on each day or in each
#' week and draws them as calendar tiles: weeks by weekday, or years by
#' week of the year.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling dates, of class \code{Date},
#'   or times of class \code{POSIXct}, counted on their date in their time
#'   zone.
#'
#' @param unit \code{Character scalar}. Counts samples per day, with weeks
#'   on the x axis and weekdays on the y axis (\code{"day"}), or per week,
#'   with the ISO week on the x axis and the year on the y axis
#'   (\code{"week"}). (Default: \code{"day"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_tile]{geom_tile}}.
#'
#' @details
#' Weeks start on Monday. Days and weeks without samples are left blank.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotCalendar
#'
#' @examples
#' # Daily sampling of 20 people; study days placed on an arbitrary
#' # calendar starting on 1 March, with no samples at weekends
#' data("temporalMicrobiome20", package = "miaTime")
#' tse <- temporalMicrobiome20
#' tse$date <- as.Date("2021-03-01") + tse$day - 1
#' tse <- tse[, !format(tse$date, "%u") %in% c("6", "7")]
#' plotCalendar(tse, time.col = "date")
#'
#' # Samples per week of the year
#' plotCalendar(tse, time.col = "date", unit = "week")
NULL

#' @rdname plotCalendar
#' @importFrom SummarizedExperiment colData
#' @importFrom stats aggregate
#' @importFrom methods setMethod
#' @export
setMethod("plotCalendar", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, unit = c("day", "week"), ...) {
        unit <- match.arg(unit)
        .require_column(x, time.col)
        date <- colData(x)[[time.col]]
        if (inherits(date, "POSIXct")) {
            date <- as.Date(date, tz = attr(date, "tzone") %||% "")
        }
        if (!inherits(date, "Date")) {
            stop("'time.col' must be a Date or POSIXct column of colData(x).",
                call. = FALSE)
        }
        date <- date[!is.na(date)]
        if (length(date) == 0L) {
            stop("'time.col' has no dates.", call. = FALSE)
        }
        # Monday of each date's week; weekday 1 = Monday.
        weekday <- as.integer(format(date, "%u"))
        monday <- date - (weekday - 1L)
        df <- if (unit == "day") {
            data.frame(.x = monday, .y = factor(weekday, levels = 7:1,
                labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")))
        } else {
            data.frame(.x = as.integer(format(monday, "%V")),
                .y = factor(format(monday, "%G")))
        }
        df <- aggregate(list(n = rep(1L, nrow(df))), by = df, FUN = sum)
        .plot_calendar(df, unit, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_tile labs
#' @importFrom utils modifyList
.plot_calendar <- function(df, unit, ...) {
    p <- ggplot(df, aes(x = .data$.x, y = .data$.y, fill = .data$n)) +
        do.call(geom_tile, modifyList(list(colour = "white"), list(...))) +
        .scale_sequential("fill") + .theme_tempoPlot() +
        labs(fill = "samples")
    if (unit == "day") {
        p + labs(x = "week starting", y = NULL)
    } else {
        p + labs(x = "week of the year", y = "year")
    }
}
