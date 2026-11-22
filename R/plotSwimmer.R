#' Plot sampling timelines of subjects
#'
#' \code{plotSwimmer()} draws one row per subject with time on the x axis
#' and one point per sample. A lane joins the samples of each subject.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   or \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object.
#'   For a \code{MultiAssayExperiment}, \code{time.col} and \code{group}
#'   are columns of its \code{colData}, and each subject row has one band
#'   per experiment that marks the samples it contains.
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
#' drawn in level order, including levels without samples. Above 60
#' subjects, subject labels are hidden.
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
#' tse <- minimalgut
#' events <- data.frame(
#'     start = c(0, 152, 432),
#'     end = c(24, NA, 436),
#'     label = c("batch", "Blautia added", "double dilution rate"))
#' plotSwimmer(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     events = events)
#'
#' # Which samples have species counts, short-chain fatty acids and
#' # optical density
#' cd <- SummarizedExperiment::colData(tse)
#' od <- !is.na(cd$OpticalDensity)
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'     experiments = list(
#'         species = tse,
#'         SCFA = SummarizedExperiment::SummarizedExperiment(
#'             t(as.matrix(cd[, c("Acetate", "Butyrate", "Propionate")]))),
#'         density = SummarizedExperiment::SummarizedExperiment(
#'             t(as.matrix(cd[od, "OpticalDensity", drop = FALSE])))),
#'     colData = cd)
#' plotSwimmer(mae, time.col = "Time.hr", group = "StudyIdentifier")
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
        .require_column(x, group)
        visit <- .check_longitudinal(x, time.col, group) == "visit"
        .check_column(x, colour.by)
        if (!visit && missing != "break" && is.null(max.gap)) {
            stop("missing = '", missing, "' requires 'max.gap' for ",
                "numeric and Date times.", call. = FALSE)
        }
        events <- .check_events(events, x, time.col, group)

        df <- as.data.frame(colData(x), optional = TRUE)
        df <- .order_time(df, time.col, group)
        if (visit && is.null(max.gap)) {
            max.gap <- 1
        }
        df <- .break_gaps(df, time.col, group, max.gap)
        missed <- .missed_visits(df, time.col, group)

        if (missing == "drop") {
            # as.character() first: c() turns a factor into its codes.
            gapped <- unique(c(as.character(df[[group]][df$.segment > 1L]),
                missed[[group]]))
            keep <- !as.character(df[[group]]) %in% gapped
            if (!any(keep)) {
                stop("All subjects have missing samples.", call. = FALSE)
            }
            if (length(gapped)) {
                message("Removed ", length(gapped),
                    " subjects with missing samples.")
            }
            df <- df[keep, , drop = FALSE]
        }

        # Rows are grouped by subject, so each run is one subject.
        runs <- rle(as.character(df[[group]]))
        first <- cumsum(runs$lengths) - runs$lengths + 1L
        key <- switch(order.by,
            first = .time_numeric(df[[time.col]])[first],
            n = -runs$lengths)
        df[[group]] <- factor(df[[group]],
            levels = rev(runs$values[order(key)]))
        .plot_swimmer(df, time.col, group, colour.by, missing,
            if (missing == "show") missed, events, ...)
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

# Dashed lines across the gaps between lane segments. Segments restart at
# 1 for each subject, so a larger segment than the previous row is a gap
# within the same subject.
.gap_lines <- function(df, time.col, group) {
    i <- which(diff(df$.segment) > 0L) + 1L
    res <- df[i, c(group, time.col), drop = FALSE]
    res$.from <- df[[time.col]][i - 1L]
    res
}

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point geom_segment
#'   labs scale_x_discrete theme element_blank
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
    # Event fills are categorical whatever colour.by is.
    p <- p + .scale_by(df, colour.by, "colour") +
        .scale_categorical("fill") + .theme_tempoPlot()
    # Labels of more than 60 rows overlap.
    if (length(rows) > 60L) {
        p <- p + theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    }
    p
}

#' @rdname plotSwimmer
#' @importFrom MultiAssayExperiment sampleMap
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotSwimmer", signature = c(x = "MultiAssayExperiment"),
    function(x, time.col, group, order.by = c("first", "n"), ...) {
        order.by <- match.arg(order.by)
        cd <- colData(x)
        for (arg in c("time.col", "group")) {
            v <- get(arg)
            if (!.is_one_name(v) || !v %in% names(cd)) {
                stop("'", arg, "' must be one column name of colData(x).",
                    call. = FALSE)
            }
        }
        time <- cd[[time.col]]
        if (!(is.numeric(time) || inherits(time, "Date") ||
                is.factor(time))) {
            stop("'time.col' must be numeric, Date or factor.",
                call. = FALSE)
        }
        map <- as.data.frame(sampleMap(x))
        exps <- levels(factor(map$assay))
        df <- data.frame(.primary = map$primary,
            .experiment = factor(map$assay, levels = exps))
        df[[time.col]] <- time[match(df$.primary, rownames(cd))]
        df[[group]] <- as.character(cd[[group]])[match(df$.primary,
            rownames(cd))]
        df <- df[!is.na(df[[group]]) & !is.na(df[[time.col]]), ,
            drop = FALSE]

        # Order subjects as the SummarizedExperiment method does.
        tnum <- .time_numeric(time)
        subject <- as.character(cd[[group]])
        first <- tapply(tnum, subject, min)
        n <- tapply(tnum, subject, length)
        key <- if (order.by == "first") first else -n
        rows <- rev(names(key)[order(key)])
        df$.row <- match(df[[group]], rows)
        k <- length(exps)
        offset <- if (k == 1L) 0 else seq(0.3, -0.3, length.out = k)
        df$.y <- df$.row + offset[as.integer(df$.experiment)]
        .plot_swimmer_mae(df, time.col, group, rows, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs
#'   scale_y_continuous scale_x_discrete
.plot_swimmer_mae <- function(df, time.col, group, rows, ...) {
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data$.y,
            colour = .data$.experiment)) +
        geom_line(aes(group = interaction(.data$.row,
            .data$.experiment)), alpha = 0.4) +
        geom_point(...) +
        scale_y_continuous(breaks = seq_along(rows), labels = rows) +
        labs(x = time.col, y = group, colour = "experiment")
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p
}
