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
#'   units. For visits, \code{NULL} breaks lanes at visits without a
#'   sample; otherwise \code{NULL} never breaks lanes. (Default: \code{NULL})
#'
#' @param missing \code{Character scalar}. How to draw missing samples:
#'   break the lane (\code{"break"}), join the samples around them with a
#'   dashed line and mark missed visits with open circles (\code{"show"}),
#'   or remove subjects with any missing sample (\code{"drop"}).
#'   (Default: \code{"break"})
#'
#' @param events \code{data.frame}. Events to draw, one per row, with
#'   columns \code{start}, and optionally \code{end}, \code{label} and the
#'   \code{group} column. \code{start} and \code{end} have the type of
#'   \code{time.col}; for visits, they are visit levels. Events with an
#'   \code{end} are drawn as bands, the others as markers. Events with a
#'   subject in the \code{group} column are drawn on the row of that
#'   subject; the others span all subjects. \code{label} sets the fill of
#'   bands and markers on a row and the line type of markers that span all
#'   subjects.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Numeric and \code{Date} times are drawn on a continuous axis. Visits are
#' drawn in level order, including levels without samples.
#'
#' For visits, a sample is missing when a subject has no sample at a visit
#' after its first one, including visits after its last sample (dropout).
#' For numeric and \code{Date} times, a gap longer than \code{max.gap} is
#' missing; \code{missing = "show"} and \code{"drop"} then require
#' \code{max.gap}.
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
#'
#' # Mark gaps longer than a year
#' plotSwimmer(tse, time.col = "time", group = "subject", max.gap = 1,
#'     missing = "show")
#'
#' # Perturbations of three bioreactors run on the same schedule
#' data("minimalgut", package = "miaTime")
#' events <- data.frame(
#'     start = c(0, 152, 432),
#'     end = c(24, NA, 436),
#'     label = c("batch", "Blautia added", "double dilution rate"))
#' plotSwimmer(minimalgut, time.col = "Time.hr", group = "StudyIdentifier",
#'     events = events)
NULL

#' @rdname plotSwimmer
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotSwimmer", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, colour.by = NULL,
            order.by = c("first", "n"), max.gap = NULL,
            missing = c("break", "show", "drop"), events = NULL, ...) {
        order.by <- match.arg(order.by)
        missing <- match.arg(missing)
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        # Errors on character and other unsupported time types.
        visit <- .time_type(x, time.col) == "visit"
        .check_column(x, colour.by)
        if (!visit && missing != "break" && is.null(max.gap)) {
            stop("missing = '", missing, "' requires 'max.gap' for ",
                "numeric and Date times.", call. = FALSE)
        }
        events <- .check_events(events, x, time.col, group)

        cols <- unique(c(time.col, group, colour.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df <- .order_time(df, time.col, group)
        if (visit && is.null(max.gap)) {
            max.gap <- 1
        }
        df <- .break_gaps(df, time.col, group, max.gap)
        missed <- .missed_visits(df, time.col, group)

        if (missing == "drop") {
            gapped <- c(df[[group]][df$.segment > 1L], missed[[group]])
            gapped <- unique(as.character(gapped))
            if (length(gapped) == length(unique(df[[group]]))) {
                stop("All subjects have missing samples.", call. = FALSE)
            }
            if (length(gapped)) {
                message("Removed ", length(gapped),
                    " subjects with missing samples.")
            }
            df <- df[!as.character(df[[group]]) %in% gapped, , drop = FALSE]
            missed <- NULL
        }

        # Rows are grouped by subject, so each run is one subject.
        runs <- rle(as.character(df[[group]]))
        first <- cumsum(runs$lengths) - runs$lengths + 1L
        key <- switch(order.by,
            first = .time_numeric(df[[time.col]])[first],
            n = -runs$lengths)
        df[[group]] <- factor(df[[group]],
            levels = rev(runs$values[order(key)]))
        if (missing != "show") {
            missed <- NULL
        }

        .plot_swimmer(df, time.col, group, colour.by, missing, missed,
            events, ...)
    }
)

# Visits without a sample after the first visit of each subject: a
# data.frame with the group and time.col columns, NULL for other times.
.missed_visits <- function(df, time.col, group) {
    time <- df[[time.col]]
    if (!is.factor(time)) {
        return(NULL)
    }
    lv <- levels(time)
    codes <- split(as.integer(time), as.character(df[[group]]))
    missed <- lapply(codes, function(cd) {
        setdiff(seq(min(cd), length(lv)), cd)
    })
    res <- data.frame(rep(names(missed), lengths(missed)),
        factor(lv[unlist(missed)], levels = lv))
    names(res) <- c(group, time.col)
    res
}

# Dashed lines across the gaps between lane segments.
.gap_lines <- function(df, time.col, group) {
    n <- nrow(df)
    i <- which(df$.segment[-1L] > df$.segment[-n] &
        df[[group]][-1L] == df[[group]][-n]) + 1L
    res <- df[i, c(group, time.col), drop = FALSE]
    res$.from <- df[[time.col]][i - 1L]
    res
}

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point geom_segment
#'   labs scale_x_discrete
#' @importFrom stats setNames
.plot_swimmer <- function(df, time.col, group, colour.by = NULL,
        missing = "break", missed = NULL, events = NULL, ...) {
    rows <- levels(df[[group]])
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data[[group]])) +
        .events_layer(events, "bands", group, rows) +
        geom_line(aes(group = interaction(.data[[group]], .data$.segment)))
    if (missing == "show") {
        gaps <- .gap_lines(df, time.col, group)
        if (nrow(gaps)) {
            p <- p + geom_segment(aes(x = .data$.from,
                xend = .data[[time.col]], yend = .data[[group]]),
                data = gaps, linetype = "dashed")
        }
    }
    p <- p +
        geom_point(aes(colour = if (!is.null(colour.by)) .data[[colour.by]]),
            ...)
    if (!is.null(missed) && nrow(missed)) {
        missed[[group]] <- factor(missed[[group]], levels = rows)
        p <- p + geom_point(data = missed, shape = 1, colour = "grey50")
    }
    p <- p + .events_layer(events, "markers", group, rows) +
        labs(x = time.col, y = group, colour = colour.by)
    # Name the legends of the event scales in use.
    mapped <- unlist(lapply(p$layers, function(l) names(l$mapping)))
    scales <- intersect(c("fill", "linetype"), mapped)
    if (length(scales)) {
        p <- p + do.call(labs, as.list(setNames(rep("event",
            length(scales)), scales)))
    }
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p
}
