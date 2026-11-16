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
#'   \code{colData(x)} that contains sampling dates, of class \code{Date}.
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
#' # calendar starting on 1 March
#' data("temporalMicrobiome20", package = "miaTime")
#' tse <- temporalMicrobiome20
#' tse$date <- as.Date("2021-03-01") + tse$day - 1
#' plotCalendar(tse, time.col = "date")
NULL

#' @rdname plotCalendar
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotCalendar", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, unit = c("day", "week"), ...) {
        unit <- match.arg(unit)
        .check_column(x, time.col)
        date <- colData(x)[[time.col]]
        if (is.null(time.col) || !inherits(date, "Date")) {
            stop("'time.col' must be a Date column of colData(x).",
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
            tab <- as.data.frame(table(date = as.character(date)),
                responseName = "n")
            d <- as.Date(as.character(tab$date))
            wd <- as.integer(format(d, "%u"))
            data.frame(.x = d - (wd - 1L), .y = factor(wd, levels = 7:1,
                labels = rev(c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat",
                    "Sun"))), n = tab$n)
        } else {
            key <- paste(format(monday, "%G"), format(monday, "%V"))
            tab <- as.data.frame(table(key), responseName = "n")
            parts <- do.call(rbind, strsplit(as.character(tab$key), " "))
            data.frame(.x = as.integer(parts[, 2L]), .y = factor(parts[, 1L]),
                n = tab$n)
        }
        .plot_calendar(df, unit, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_tile labs scale_fill_gradient
#'   scale_x_date
.plot_calendar <- function(df, unit, ...) {
    p <- ggplot(df, aes(x = .data$.x, y = .data$.y, fill = .data$n)) +
        geom_tile(colour = "white", ...) +
        scale_fill_gradient(low = "#c6dbef", high = "#08519c") +
        labs(fill = "samples")
    if (unit == "day") {
        p + labs(x = "week starting", y = NULL)
    } else {
        p + labs(x = "week of the year", y = "year")
    }
}
