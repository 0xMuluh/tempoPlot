#' Plot longitudinal trends
#'
#' \code{plotTrend()} plots a feature, colData column or ordination
#' component over time: one line per subject or boxes per time point, with
#' optional group summaries, model fits, and values relative to each
#' subject's first sample or to an event.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param y \code{Character scalar}. Specifies the feature, column or
#'   ordination component to plot: a column of \code{colData(x)}, a row of
#'   \code{assay(x, assay.type)}, or \code{"<dimred>:<component>"}.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits. Time differences, as in
#'   \code{max.gap} and \code{period}, are in days for \code{Date} and
#'   hours for \code{POSIXct} times.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to
#'   retrieve \code{y} from when \code{y} is a feature name.
#'   (Default: \code{"counts"})
#'
#' @param layout \code{Character scalar}. Specifies the plot layout:
#'   \code{"line"} draws one line per subject; \code{"box"} and
#'   \code{"violin"} draw boxes or violins of \code{y} per time point, one
#'   per \code{group.by} group side by side. Subject lines and points stay
#'   at the centre of each time point. (Default: \code{"line"})
#'
#' @param show.individuals \code{Logical scalar}. When \code{TRUE}, draws
#'   subject lines and points. When \code{FALSE}, hides subject lines and
#'   points, leaving only boxes or the smooth curve. (Default: \code{TRUE})
#'
#' @param smooth \code{Character scalar}. Specifies a group summary on top:
#'   \code{"none"}, \code{"mean"} or \code{"median"} (over the samples at
#'   each distinct time), or \code{"loess"} (loess curve). With irregular
#'   times, where most times have one sample, use \code{"loess"} or a
#'   factor of visits. (Default: \code{"none"})
#'
#' @param band \code{Numeric vector} of length 2. Draws the range between
#'   these quantiles of \code{y} at each distinct time, such as
#'   \code{c(0.25, 0.75)} for the interquartile range, per
#'   \code{group.by} group. (Default: \code{NULL})
#'
#' @param fit Model results on the scale of \code{y}, drawn over the data:
#'   a \code{character scalar} naming a numeric \code{colData(x)} column of
#'   fitted values per sample, drawn as a dashed line per subject, or a
#'   \code{data.frame} with the columns \code{time.col} and
#'   \code{estimate}, optionally \code{lower} and \code{upper} for an
#'   interval and the \code{group.by} column for one curve per group.
#'   Without the \code{group.by} column, the table draws one curve, such
#'   as a microbiota-for-age reference. Not combined with \code{baseline}
#'   or \code{rolling}. For effects, such as differences between groups,
#'   use \code{\link{plotTimeEffect}}. (Default: \code{NULL})
#'
#' @param baseline \code{Character scalar}. Expresses \code{y} relative to
#'   each subject's first sample: unchanged (\code{"none"}), as the
#'   difference (\code{"difference"}), the ratio (\code{"fold"}) or its log2
#'   (\code{"log2fold"}). Ratios need positive first values.
#'   (Default: \code{"none"})
#'
#' @param rolling \code{Character scalar}. Replaces \code{y} by a rolling
#'   statistic over the last \code{window} samples of each subject:
#'   \code{"mean"}, \code{"variance"}, \code{"acf1"} (lag-1
#'   autocorrelation) or \code{"slope"} (change of \code{y} per
#'   \code{time.col} unit, from a least-squares line; per visit step for
#'   visits). Rising variance and autocorrelation are early warning signals
#'   of a transition. Windows count samples, also across breaks of the
#'   line; subjects with fewer than \code{window} samples are dropped.
#'   (Default: \code{"none"})
#'
#' @param window \code{Integer scalar}. Number of samples in the rolling
#'   window. (Default: \code{5})
#'
#' @param detection \code{Numeric scalar}. Splits \code{y} into two panels,
#'   as in two-part models of zero-inflated data: its presence, a value
#'   above \code{detection}, as 1 or 0 (\code{"Prevalence"}), and its values
#'   in the samples where it is present (\code{"When present"}). Subject
#'   lines of the values break at absences, and \code{smooth = "mean"} gives
#'   the prevalence at each time. Not combined with \code{baseline},
#'   \code{rolling}, \code{fit} or boxes and violins. (Default: \code{NULL})
#'
#' @param period \code{Numeric scalar}. Wraps time into cycles of this
#'   length, in \code{time.col} units, to show daily or seasonal patterns;
#'   each cycle is drawn as its own line. Cycles start at time 0, which for
#'   \code{Date} and \code{POSIXct} times is 1970-01-01 (UTC), a Thursday.
#'   Not combined with \code{events} or a \code{fit} table.
#'   (Default: \code{NULL})
#'
#' @param events \code{data.frame}. Events or intervals to mark, one per
#'   row, with columns \code{start} and optionally \code{end} and
#'   \code{label}, in the units of \code{time.col} (visit levels for
#'   visits), relative to the event with \code{align.to}: interventions
#'   such as antibiotic courses, or the time intervals where groups differ,
#'   as estimated by MetaLonDA or SplinectomeR. Rows with an \code{end} are
#'   drawn as grey bands with their label at the top, the others as
#'   vertical lines. Events apply to all subjects, unless they have the
#'   \code{group} column and \code{facet.by} is \code{group}: then each
#'   event is drawn in the panel of its subject. To draw time relative to an
#'   event of each subject, use \code{align.to}. (Default: \code{NULL})
#'
#' @param align.to \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} with the time of an event for each subject, such as
#'   the start of a treatment, in the units of \code{time.col}. Time is then
#'   drawn relative to the event, marked by a dashed line at 0, and
#'   \code{smooth}, \code{band} and \code{fit} are drawn separately before
#'   the event and from time 0 on. The column must not contain NA values, so
#'   remove subjects without the event first. Not available for visits.
#'   (Default: \code{NULL})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that colours subjects, smooth curves and fits and
#'   splits boxes and violins. (Default: \code{NULL})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} to facet by. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Breaks a subject's line where
#'   consecutive samples are further apart than this, in \code{time.col}
#'   units. For visits, \code{NULL} breaks lines at visits without a
#'   sample; otherwise \code{NULL} never breaks lines. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to the main layer
#'   (\code{\link[ggplot2:geom_line]{geom_line}} for \code{layout = "line"},
#'   \code{\link[ggplot2:geom_boxplot]{geom_boxplot}} for
#'   \code{layout = "box"} or \code{\link[ggplot2:geom_violin]{geom_violin}}
#'   for \code{layout = "violin"}).
#'
#' @details
#' Numeric and \code{Date} times are drawn on a continuous axis. Visits are
#' drawn in level order, including levels without samples. With
#' \code{align.to}, the times of \code{fit} tables and \code{events} are
#' relative to the event as well; for \code{Date} times they are in days.
#' The plot data hold all columns of \code{colData(x)}, so further layers
#' can use them.
#'
#' Fits that stop at an aligned event follow DiGiulio et al. (2015,
#' doi:10.1073/pnas.1502875112, Fig. 5). The prevalence and non-zero panels
#' of \code{detection} follow MaAsLin 3 (Nickols et al. 2026,
#' doi:10.1038/s41592-025-02923-9, Fig. 1b). Cycles of varying length, such
#' as menstrual cycles, are compared by rescaling each cycle to a common
#' length in a time column, then wrapping it with \code{period} (Gajer et
#' al. 2012, doi:10.1126/scitranslmed.3003605).
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Bokulich NA, et al. (2018). q2-longitudinal: longitudinal and
#' paired-sample analyses of microbiome data. \emph{mSystems} 3(6):e00219-18.
#'
#' Dakos V, et al. (2012). Methods for detecting early warnings of critical
#' transitions in time series illustrated using simulated ecological data.
#' \emph{PLoS One} 7(7):e41010.
#'
#' DiGiulio DB, et al. (2015). Temporal and spatial variation of the human
#' microbiota during pregnancy. \emph{PNAS} 112(35):11060-11065.
#'
#' Gajer P, et al. (2012). Temporal dynamics of the human vaginal
#' microbiota. \emph{Sci Transl Med} 4(132):132ra52.
#'
#' Nickols WA, et al. (2026). MaAsLin 3: refining and extending generalized
#' multivariable linear models for meta-omic association discovery.
#' \emph{Nat Methods} 23(3):554-564.
#'
#' Subramanian S, et al. (2014). Persistent gut microbiota immaturity in
#' malnourished Bangladeshi children. \emph{Nature} 510(7505):417-421.
#'
#' @seealso \code{\link{plotTimeEffect}} for effects over time.
#'
#' @name plotTrend
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' counts <- SummarizedExperiment::assay(tse, "counts")
#' tse <- tse[, !duplicated(t(counts))]
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
#'
#' # Divergence from each subject's first sample, computed by miaTime,
#' # per visit (the order of each subject's samples), by sex
#' tse <- miaTime::addBaselineDivergence(tse, group = "subject",
#'     time.col = "time", assay.type = "counts")
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' plotTrend(tse, y = "divergence", time.col = "visit", group = "subject",
#'     group.by = "sex", layout = "box")
#'
#' # Divergence from each subject's previous sample, computed by miaTime
#' tse <- miaTime::addStepwiseDivergence(tse, group = "subject",
#'     time.col = "time", assay.type = "counts",
#'     name = c("stepwise", "step_time", "step_ref"))
#' plotTrend(tse, y = "stepwise", time.col = "visit", group = "subject",
#'     group.by = "sex", layout = "box")
#'
#' # Change in diversity from the first visit, mean by sex
#' plotTrend(tse, y = "diversity", time.col = "visit", group = "subject",
#'     group.by = "sex", baseline = "difference", smooth = "mean")
#'
#' # Median diversity with the interquartile range, by sex
#' plotTrend(tse, y = "diversity", time.col = "visit", group = "subject",
#'     group.by = "sex", smooth = "median", band = c(0.25, 0.75),
#'     show.individuals = FALSE)
#'
#' # Mean diversity per visit and sex from a linear model, with 95%
#' # confidence intervals; each subject has one sample per visit
#' cd <- as.data.frame(SummarizedExperiment::colData(tse))
#' model <- lm(diversity ~ visit * sex, data = cd)
#' means <- expand.grid(visit = levels(tse$visit), sex = levels(tse$sex))
#' means <- cbind(means, predict(model, means, interval = "confidence"))
#' names(means)[3:5] <- c("estimate", "lower", "upper")
#' plotTrend(tse, y = "diversity", time.col = "visit", group = "subject",
#'     group.by = "sex", fit = means)
#'
#' # Acetate in three replicate bioreactors sampled every few hours
#' data("minimalgut", package = "miaTime")
#' tse <- minimalgut
#' plotTrend(tse, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier")
#'
#' # Fitted acetate of each bioreactor from a cubic time trend per reactor
#' cd <- as.data.frame(SummarizedExperiment::colData(tse))
#' model <- lm(Acetate ~ poly(Time_hr_num, 3) * StudyIdentifier, data = cd)
#' tse$fitted <- fitted(model)
#' plotTrend(tse, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier",
#'     fit = "fitted")
#'
#' # Hours since Blautia hydrogenotrophica was added at 152 h, against the
#' # range of acetate before it was added, with loess curves fitted
#' # separately before and after
#' tse$added <- 152
#' before <- range(tse$Acetate[tse$Time_hr_num < 152])
#' plotTrend(tse, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier",
#'     align.to = "added", smooth = "loess") +
#'     ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
#'         ymin = before[1L], ymax = before[2L], alpha = 0.2)
#'
#' # Share of bioreactors where Blautia was detected, and its counts when
#' # detected
#' plotTrend(tse, y = "Blautia_hydrogenotrophica", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", align.to = "added", detection = 0,
#'     smooth = "mean")
#'
#' # Acetate over the day, one line per day and bioreactor
#' plotTrend(tse, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier", period = 24)
#'
#' # Rolling variance of acetate over 10 samples
#' plotTrend(tse, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier",
#'     rolling = "variance", window = 10)
#'
#' # Acetate with the batch phase and the addition of Blautia marked
#' events <- data.frame(start = c(0, 152), end = c(24, NA),
#'     label = c("batch", "Blautia added"))
#' plotTrend(tse, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier",
#'     events = events)
NULL

#' @rdname plotTrend
#' @importFrom ggplot2 labs
#' @importFrom methods setMethod
#' @export
setMethod("plotTrend", signature = c(x = "SummarizedExperiment"),
    function(x, y, time.col, group = NULL, assay.type = "counts",
            layout = c("line", "box", "violin"), show.individuals = TRUE,
            smooth = c("none", "mean", "median", "loess"), band = NULL,
            fit = NULL,
            baseline = c("none", "difference", "fold", "log2fold"),
            rolling = c("none", "mean", "variance", "acf1", "slope"),
            window = 5, detection = NULL, period = NULL, events = NULL,
            align.to = NULL,
            group.by = NULL, facet.by = NULL, max.gap = NULL, ...) {
        layout <- match.arg(layout)
        smooth <- match.arg(smooth)
        baseline <- match.arg(baseline)
        rolling <- match.arg(rolling)
        type <- .check_longitudinal(x, time.col, group)
        .check_column(x, group.by)
        .check_column(x, facet.by)
        x <- .align_time(x, time.col, group, align.to)
        events <- .check_events(events, x, time.col)
        .check_trend_args(x, time.col, group, type, layout,
            show.individuals, band, fit, baseline, rolling, window,
            detection, period, events, facet.by)
        df <- .trend_data(x, y, time.col, group, assay.type, baseline,
            rolling, window, period, max.gap)
        if (!is.null(detection)) {
            df <- .two_part(df, detection)
        }

        ylab <- switch(baseline, none = y,
            difference = paste("Change in", y),
            fold = paste("Fold change in", y),
            log2fold = paste("log2 fold change in", y))
        if (rolling != "none") {
            ylab <- paste0("Rolling ", rolling, " of ", ylab, " (", window,
                " samples)")
        }
        xlab <- paste0(time.col, if (!is.null(align.to))
            paste(" since", align.to), if (!is.null(period))
            paste(" modulo", period))
        .plot_trend(df, time.col, layout, show.individuals, smooth, band,
            fit, group.by, facet.by, events, !is.null(align.to), ...) +
            labs(x = xlab, y = ylab)
    }
)

# Errors on the arguments of plotTrend that its helpers do not check, and
# on the combinations it cannot draw.
#' @importFrom SummarizedExperiment colData
.check_trend_args <- function(x, time.col, group, type, layout,
        show.individuals, band, fit, baseline, rolling, window, detection,
        period, events, facet.by) {
    .check_flag(show.individuals)
    if (!is.null(band) && !(is.numeric(band) && length(band) == 2L &&
            all(band >= 0 & band <= 1) && band[1L] < band[2L])) {
        stop("'band' must be NULL or two increasing probabilities.",
            call. = FALSE)
    }
    if (!(is.null(fit) ||
            .is_one_name(fit) && is.numeric(colData(x)[[fit]]) ||
            is.data.frame(fit) &&
                all(c(time.col, "estimate") %in% names(fit)))) {
        stop("'fit' must be a numeric column of colData(x) or a data.frame ",
            "with the columns '", time.col, "' and estimate.", call. = FALSE)
    }
    if (!is.null(fit) && (baseline != "none" || rolling != "none")) {
        stop("'fit' cannot be combined with 'baseline' or 'rolling'.",
            call. = FALSE)
    }
    .check_number(window, 2)
    .check_number(detection, null = TRUE)
    if (!is.null(detection) && (layout != "line" || baseline != "none" ||
            rolling != "none" || !is.null(fit))) {
        stop("'detection' needs layout = \"line\", without 'baseline', ",
            "'rolling' or 'fit'.", call. = FALSE)
    }
    .check_number(period, 0, open = TRUE, null = TRUE)
    if (!is.null(period) && (type == "visit" || !is.null(events) ||
            is.data.frame(fit))) {
        stop("'period' needs numeric, Date or POSIXct times, without ",
            "'events' or a 'fit' table.", call. = FALSE)
    }
    if (any(group %in% names(events)) && !identical(facet.by, group)) {
        stop("'events' of single subjects need facet.by = group; use ",
            "'align.to' for an event time per subject.", call. = FALSE)
    }
}

# x with time.col made relative to the event time in the align.to column,
# which must have the type of time.col and be constant within subjects.
#' @importFrom SummarizedExperiment colData
.align_time <- function(x, time.col, group, align.to) {
    .check_column(x, align.to)
    if (is.null(align.to)) {
        return(x)
    }
    time <- colData(x)[[time.col]]
    event <- colData(x)[[align.to]]
    if (is.factor(time) || !(if (is.numeric(time)) is.numeric(event) else
            inherits(event, class(time)[1L]))) {
        stop("'align.to' must be a column with the type of 'time.col', ",
            "which must be numeric, Date or POSIXct.", call. = FALSE)
    }
    if (anyNA(event)) {
        stop("'align.to' must not contain NA values; remove the ",
            "subjects without the event.", call. = FALSE)
    }
    subject <- if (is.null(group)) "" else colData(x)[[group]]
    if (anyDuplicated(unique(data.frame(subject, event))$subject)) {
        stop("'align.to' must be constant within each subject.",
            call. = FALSE)
    }
    x[[time.col]] <- .time_numeric(time) - .time_numeric(event)
    x
}

# colData of x ordered by time, with the value .y (relative to the first
# sample with baseline, or a rolling statistic), the subject .id and one
# .line per subject, gap segment and period cycle.
#' @importFrom SummarizedExperiment colData
#' @importFrom stats ave
.trend_data <- function(x, y, time.col, group, assay.type, baseline,
        rolling, window, period, max.gap) {
    df <- as.data.frame(colData(x), optional = TRUE)
    df$.y <- unname(.resolve_y(x, y, assay.type))
    df <- .break_gaps(.order_time(df, time.col, group), time.col, group,
        max.gap)
    df$.id <- if (is.null(group)) 1L else df[[group]]
    if (baseline != "none") {
        first <- ave(df$.y, df$.id, FUN = function(v) v[1L])
        if (baseline != "difference" && any(first <= 0, na.rm = TRUE)) {
            stop("'baseline' ratios need positive first values; use ",
                "\"difference\".", call. = FALSE)
        }
        df$.y <- switch(baseline, difference = df$.y - first,
            fold = df$.y / first, log2fold = log2(df$.y / first))
    }
    if (rolling != "none") {
        tnum <- .time_numeric(df[[time.col]])
        v <- df$.y
        for (i in split(seq_along(v), df$.id)) {
            df$.y[i] <- .rolling(v[i], rolling, as.integer(window), tnum[i])
        }
        # The first window - 1 samples of each subject have no value.
        df <- df[!is.na(df$.y), , drop = FALSE]
        if (nrow(df) == 0L) {
            stop("No subject has 'window' samples for 'rolling'.",
                call. = FALSE)
        }
    }
    cycle <- 0L
    if (!is.null(period)) {
        t <- .time_numeric(df[[time.col]])
        cycle <- t %/% period
        df[[time.col]] <- t %% period
    }
    df$.line <- interaction(df$.id, df$.segment, cycle, drop = TRUE)
    df
}

# df twice, in the panels .part: presence of .y above detection as 1 or 0,
# and .y of the samples where it is present, with lines broken at absences.
.two_part <- function(df, detection) {
    presence <- df
    presence$.y <- as.numeric(df$.y > detection)
    presence$.part <- "Prevalence"
    present <- presence$.y %in% 1
    df$.line <- interaction(df$.line, cumsum(!present), drop = TRUE)
    df$.part <- "When present"
    df <- rbind(presence, df[present, , drop = FALSE])
    df$.part <- factor(df$.part, c("Prevalence", "When present"))
    df
}

#' @importFrom utils modifyList
#' @importFrom ggplot2 ggplot aes .data geom_line geom_point geom_boxplot
#'   geom_violin geom_vline facet_wrap facet_grid vars labs
#'   scale_x_discrete
.plot_trend <- function(df, time.col, layout, show.individuals, smooth,
        band, fit, group.by, facet.by, events, aligned, ...) {
    # After align.to, summaries and fits stop at the event: one curve
    # before it and one from time 0 on.
    df$.side <- if (aligned) df[[time.col]] >= 0 else TRUE
    # group.by sets the fill of boxes, violins and bands too.
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data$.y,
            colour = if (!is.null(group.by)) .data[[group.by]],
            fill = if (!is.null(group.by)) .data[[group.by]])) +
        .events_layer(events, "bands", text = TRUE, x = df[[time.col]])
    dots <- list(...)
    if (layout != "line") {
        # One box or violin per time and group.by group; individual points
        # replace the outlier points of boxes.
        geom <- if (layout == "box") geom_boxplot else geom_violin
        args <- if (layout == "box") list(outliers = !show.individuals)
        p <- p + do.call(geom, c(list(aes(group = if (is.null(group.by))
                .data[[time.col]] else
                interaction(.data[[time.col]], .data[[group.by]]),
            colour = NULL)), modifyList(as.list(args), dots)))
    }
    # Individuals fade under a summary, and when there are many of them.
    summarised <- layout != "line" || smooth != "none" || !is.null(fit) ||
        !is.null(band)
    alpha <- if (summarised) 0.3 else min(1, 30 / nlevels(df$.line))
    fade <- if (alpha < 1) list(alpha = alpha) else list()
    if (show.individuals) {
        line_args <- if (layout == "line") modifyList(fade, dots) else fade
        p <- p + do.call(geom_line, c(list(aes(group = .data$.line)),
            line_args)) + do.call(geom_point, fade)
    }
    p + list(.trend_summaries(df, time.col, smooth, band, group.by),
        if (is.character(fit)) geom_line(aes(y = .data[[fit]],
            group = interaction(.data$.line, .data$.side)),
            linetype = "dashed", linewidth = 1),
        if (is.data.frame(fit)) .fit_layers(fit, time.col, group.by,
            split = aligned),
        .events_layer(events, "markers"),
        if (aligned) geom_vline(xintercept = 0, linetype = "dashed",
            colour = "grey40"),
        if (!is.null(df[[".part"]])) facet_grid(rows = vars(.data$.part),
            cols = if (!is.null(facet.by)) vars(.data[[facet.by]]),
            scales = "free_y") else if (!is.null(facet.by))
            facet_wrap(vars(.data[[facet.by]])),
        if (is.factor(df[[time.col]])) scale_x_discrete(drop = FALSE),
        labs(colour = group.by, fill = group.by),
        .scale_by(df, group.by, c("colour", "fill")),
        .theme_tempoPlot())
}

# Group summaries over the samples at each time: the band between two
# quantiles, and the mean, median or loess curve.
#' @importFrom ggplot2 aes .data stat_summary geom_smooth
#' @importFrom stats median quantile
.trend_summaries <- function(df, time.col, smooth, band, group.by) {
    if ((smooth %in% c("mean", "median") || !is.null(band)) &&
            !is.factor(df[[time.col]]) &&
            median(table(paste(df[[time.col]], df[[".part"]]))) < 2) {
        warning("Most times have a single sample, so 'smooth' and 'band' ",
            "summarise single samples; use smooth = \"loess\" or a factor ",
            "of visits.", call. = FALSE)
    }
    by <- aes(group = interaction(if (is.null(group.by)) 1L else
        .data[[group.by]], .data$.side))
    list(if (!is.null(band)) stat_summary(by,
            fun.min = function(v) quantile(v, band[1L], names = FALSE),
            fun.max = function(v) quantile(v, band[2L], names = FALSE),
            fun = median, geom = "ribbon", alpha = 0.2, colour = NA),
        switch(smooth, none = NULL,
            loess = geom_smooth(by, method = "loess", formula = y ~ x,
                se = FALSE, linewidth = 1),
            stat_summary(by, fun = if (smooth == "mean") mean else median,
                geom = "line", linewidth = 1)))
}

# Line and interval of a results table, one per group.by group if the
# table has that column, split at time 0 with split. ... goes to the line.
# Shared with plotTimeEffect.
#' @importFrom utils modifyList
#' @importFrom ggplot2 aes .data geom_line geom_ribbon labs
.fit_layers <- function(fit, time.col, group.by, split = FALSE, ...) {
    by <- if (!is.null(group.by) && group.by %in% names(fit)) group.by
    fit$.curve <- interaction(if (is.null(by)) 1L else fit[[by]],
        if (split) fit[[time.col]] >= 0 else TRUE)
    line <- do.call(geom_line, c(list(aes(x = .data[[time.col]],
            y = .data$estimate, group = .data$.curve,
            colour = if (!is.null(by)) .data[[by]]),
        data = fit, inherit.aes = FALSE),
        modifyList(list(linewidth = 1), list(...))))
    if (!all(c("lower", "upper") %in% names(fit))) {
        return(list(line))
    }
    list(geom_ribbon(aes(x = .data[[time.col]], ymin = .data$lower,
            ymax = .data$upper, group = .data$.curve,
            fill = if (!is.null(by)) .data[[by]]),
        data = fit, inherit.aes = FALSE, alpha = 0.2), line, labs(fill = by))
}

# Statistic over the last window values; NA until the window is full.
# t gives the times, used by "slope".
#' @importFrom stats var cor cov
.rolling <- function(v, stat, window, t) {
    n <- length(v)
    out <- rep(NA_real_, n)
    for (i in seq_len(n)[seq_len(n) >= window]) {
        j <- (i - window + 1L):i
        w <- v[j]
        # A constant window has no autocorrelation: cor() warns and gives
        # NA, which drops the point like an incomplete window.
        out[i] <- switch(stat, mean = mean(w, na.rm = TRUE),
            variance = var(w, na.rm = TRUE),
            acf1 = suppressWarnings(cor(w[-1L], w[-window],
                use = "complete.obs")),
            slope = cov(t[j], w, use = "complete.obs") /
                var(t[j], na.rm = TRUE))
    }
    out
}
