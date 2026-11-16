#' Plot longitudinal trends
#'
#' \code{plotTrend()} plots the trajectory of a feature, colData column or
#' ordination component over time for each subject, with optional group
#' summaries, boxplot paired layouts, and values relative to each
#' subject's first sample.
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
#'   \code{"line"} draws one line per subject; \code{"box"} and
#'   \code{"violin"} draw boxes or violins of \code{y} per time point;
#'   \code{"ridge"} draws the density of \code{y} at each time point as
#'   stacked ridges (needs ggridges).
#'   (Default: \code{"line"})
#'
#' @param show.individuals \code{Logical scalar}. When \code{TRUE}, draws
#'   subject lines and points. When \code{FALSE}, hides subject lines and
#'   points, leaving only boxes or the smooth curve. (Default: \code{TRUE})
#'
#' @param smooth \code{Character scalar}. Specifies a group summary on top:
#'   \code{"none"}, \code{"mean"} or \code{"median"} (at each time point),
#'   or \code{"loess"} (loess curve). (Default: \code{"none"})
#'
#' @param band \code{Numeric vector} of length 2. Draws the range between
#'   these quantiles of \code{y} at each time point, such as
#'   \code{c(0.25, 0.75)} for the interquartile range, per
#'   \code{group.by} group. (Default: \code{NULL})
#'
#' @param fit Model results to draw over the data: the name of a numeric
#'   \code{colData(x)} column of fitted values per sample, drawn as a dashed
#'   line per subject, or a table with the columns \code{time.col} and
#'   \code{estimate}, optionally \code{lower} and \code{upper} for an
#'   interval and the \code{group.by} column for groups. Not combined with
#'   \code{baseline}. (Default: \code{NULL})
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
#'   of a transition. (Default: \code{"none"})
#'
#' @param window \code{Integer scalar}. Number of samples in the rolling
#'   window. (Default: \code{5})
#'
#' @param period \code{Numeric scalar}. Wraps time into cycles of this
#'   length, in \code{time.col} units (days for \code{Date} times), to
#'   show daily or seasonal patterns; each cycle is drawn as its own line.
#'   (Default: \code{NULL})
#'
#' @param reference A reference curve to draw behind the data, such as
#'   microbiota-for-age: a table with the columns \code{time.col} and
#'   \code{estimate}, optionally \code{lower} and \code{upper} for a band.
#'   (Default: \code{NULL})
#'
#' @param limits \code{Numeric vector} of length 2. Draws a band between
#'   these values of \code{y}, such as a normal range.
#'   (Default: \code{NULL})
#'
#' @param align.to \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} with the time of an event for each subject, such as
#'   the start of a treatment, in the units of \code{time.col}. Time is then
#'   drawn relative to the event. Not available for visits.
#'   (Default: \code{NULL})
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
#'   (\code{\link[ggplot2:geom_line]{geom_line}} for \code{layout = "line"},
#'   \code{\link[ggplot2:geom_violin]{geom_violin}} for
#'   \code{layout = "violin"}
#'   or \code{\link[ggplot2:geom_boxplot]{geom_boxplot}} for
#'   \code{layout = "box"}).
#'
#' @details
#' Numeric and \code{Date} times are drawn on a continuous axis. Visits are
#' drawn in level order, including levels without samples. With
#' \code{align.to}, the times of \code{fit} and \code{reference} tables are
#' relative to the event as well; for \code{Date} times they are in days.
#'
#' @return
#' A \code{ggplot} object.
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
#' # per visit (the order of each subject's samples)
#' tse <- miaTime::addBaselineDivergence(tse, group = "subject",
#'     time.col = "time", assay.type = "counts")
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' plotTrend(tse, y = "divergence", time.col = "visit", group = "subject",
#'     layout = "box")
#'
#' # Change in diversity from the first visit, mean by sex
#' plotTrend(tse, y = "diversity", time.col = "visit", group = "subject",
#'     group.by = "sex", baseline = "difference", smooth = "mean")
#'
#' # Mean diversity per visit and sex from a linear model, with 95%
#' # confidence intervals
#' cd <- as.data.frame(SummarizedExperiment::colData(tse))
#' model <- lm(diversity ~ visit * sex, data = cd)
#' means <- expand.grid(visit = levels(tse$visit), sex = levels(tse$sex))
#' means <- cbind(means, predict(model, means, interval = "confidence"))
#' names(means)[3:5] <- c("estimate", "lower", "upper")
#' plotTrend(tse, y = "diversity", time.col = "visit", group = "subject",
#'     group.by = "sex", fit = means)
#'
#' # Fitted acetate of each bioreactor from a cubic time trend per reactor
#' data("minimalgut", package = "miaTime")
#' cd <- as.data.frame(SummarizedExperiment::colData(minimalgut))
#' model <- lm(Acetate ~ poly(Time_hr_num, 3) * StudyIdentifier, data = cd)
#' minimalgut$fitted <- fitted(model)
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier",
#'     fit = "fitted")
#'
#' # Acetate in three replicate bioreactors sampled every few hours
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier")
#'
#' # Hours since Blautia hydrogenotrophica was added at 152 h, against the
#' # range of acetate before it was added
#' minimalgut$added <- 152
#' before <- minimalgut$Time_hr_num < 152
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier",
#'     align.to = "added", limits = range(minimalgut$Acetate[before]))
#'
#' # Acetate over the day, one line per day and bioreactor
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier", period = 24)
#'
#' # Rolling variance of acetate over 10 samples
#' plotTrend(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier", group.by = "StudyIdentifier",
#'     rolling = "variance", window = 10)
#'
#' # Median diversity with the interquartile range, by sex
#' plotTrend(tse, y = "diversity", time.col = "visit", group = "subject",
#'     group.by = "sex", smooth = "median", band = c(0.25, 0.75),
#'     show.individuals = FALSE)
#'
#' # Distribution of diversity at each visit
#' plotTrend(tse, y = "diversity", time.col = "visit", group = "subject",
#'     layout = "ridge")
NULL

#' @rdname plotTrend
#' @importFrom SummarizedExperiment colData
#' @importFrom stats ave
#' @importFrom methods setMethod
#' @export
setMethod("plotTrend", signature = c(x = "SummarizedExperiment"),
    function(x, y, time.col, group = NULL, assay.type = NULL,
            layout = c("line", "box", "violin", "ridge"),
            show.individuals = TRUE,
            smooth = c("none", "mean", "median", "loess"), band = NULL,
            fit = NULL,
            baseline = c("none", "difference", "fold", "log2fold"),
            rolling = c("none", "mean", "variance", "acf1", "slope"),
            window = 5,
            period = NULL, reference = NULL, limits = NULL, align.to = NULL,
            group.by = NULL, facet.by = NULL, max.gap = NULL, ...) {
        layout <- match.arg(layout)
        smooth <- match.arg(smooth)
        baseline <- match.arg(baseline)
        rolling <- match.arg(rolling)
        if (!(is.numeric(window) && length(window) == 1L && window >= 2)) {
            stop("'window' must be a single number of at least 2.",
                call. = FALSE)
        }
        if (!is.null(period) && !(is.numeric(period) &&
                length(period) == 1L && period > 0)) {
            stop("'period' must be NULL or a single positive number.",
                call. = FALSE)
        }
        if (!is.null(band) && !(is.numeric(band) && length(band) == 2L &&
                all(band >= 0 & band <= 1) && band[1L] < band[2L])) {
            stop("'band' must be NULL or two increasing probabilities.",
                call. = FALSE)
        }
        if (layout == "ridge" &&
                !requireNamespace("ggridges", quietly = TRUE)) {
            stop("layout = \"ridge\" requires the ggridges package.",
                call. = FALSE)
        }
        if (!(isTRUE(show.individuals) || isFALSE(show.individuals))) {
            stop("'show.individuals' must be TRUE or FALSE.", call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        # Errors on character and other unsupported time types.
        type <- .time_type(x, time.col)
        if (!is.null(period) && type == "visit") {
            stop("'period' needs numeric or Date times.", call. = FALSE)
        }
        .check_column(x, group.by)
        .check_column(x, facet.by)
        .check_column(x, align.to)
        if (!is.null(align.to)) {
            event <- colData(x)[[align.to]]
            if (type == "visit" || !(is.numeric(event) && type == "numeric" ||
                    inherits(event, "Date") && type == "date")) {
                stop("'align.to' must be a column with the type of ",
                    "'time.col', which must be numeric or Date.",
                    call. = FALSE)
            }
        }
        if (!is.null(reference)) {
            reference <- as.data.frame(reference)
            if (!all(c(time.col, "estimate") %in% names(reference))) {
                stop("'reference' must be a table with the columns '",
                    time.col, "' and estimate.", call. = FALSE)
            }
        }
        if (!is.null(limits) && !(is.numeric(limits) &&
                length(limits) == 2L && !anyNA(limits))) {
            stop("'limits' must be two numbers.", call. = FALSE)
        }
        if (is.character(fit)) {
            .check_column(x, fit)
            if (!is.numeric(colData(x)[[fit]])) {
                stop("'fit' must be a numeric column of colData(x).",
                    call. = FALSE)
            }
        } else if (!is.null(fit)) {
            fit <- as.data.frame(fit)
            if (!all(c(time.col, "estimate") %in% names(fit))) {
                stop("'fit' must be a colData column or a table with the ",
                    "columns '", time.col, "' and estimate.", call. = FALSE)
            }
        }
        if (!is.null(fit) && baseline != "none") {
            stop("'fit' cannot be combined with 'baseline'.", call. = FALSE)
        }

        cols <- unique(c(time.col, group, group.by, facet.by, align.to,
            if (is.character(fit)) fit))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        if (!is.null(align.to)) {
            ids <- if (is.null(group)) rep(1L, nrow(df)) else df[[group]]
            event <- as.numeric(df[[align.to]])
            if (any(tapply(event, ids, function(e)
                    length(unique(e)) > 1L))) {
                stop("'align.to' must be constant within each subject.",
                    call. = FALSE)
            }
            df[[time.col]] <- as.numeric(df[[time.col]]) - event
        }
        df$.y <- unname(.resolve_y(x, y, assay.type))
        if (is.character(fit)) {
            df$.fit <- df[[fit]]
        }
        df <- .order_time(df, time.col, group)
        df <- .break_gaps(df, time.col, group, max.gap)

        ids <- if (is.null(group)) rep(1L, nrow(df)) else df[[group]]
        if (baseline != "none") {
            first <- ave(df$.y, ids, FUN = function(v) v[1L])
            if (baseline != "difference" && any(first <= 0, na.rm = TRUE)) {
                stop("'baseline' ratios need positive first values; use ",
                    "\"difference\".", call. = FALSE)
            }
            df$.y <- switch(baseline, difference = df$.y - first,
                fold = df$.y / first, log2fold = log2(df$.y / first))
        }
        if (rolling != "none") {
            tnum <- .time_numeric(df[[time.col]])
            y <- df$.y
            df$.y <- ave(as.numeric(seq_along(y)), ids, FUN = function(i) {
                i <- as.integer(i)
                .rolling(y[i], rolling, as.integer(window), tnum[i])
            })
            # The first window - 1 samples of each subject have no value.
            keep <- !is.na(df$.y)
            df <- df[keep, , drop = FALSE]
            ids <- ids[keep]
        }
        cycle <- 0L
        if (!is.null(period)) {
            t <- as.numeric(df[[time.col]])
            cycle <- floor(t / period)
            df[[time.col]] <- t %% period
        }
        df$.line <- interaction(ids, df$.segment, cycle, drop = TRUE)

        xlab <- if (is.null(align.to)) time.col else
            paste(time.col, "since", align.to)
        if (!is.null(period)) {
            xlab <- paste(xlab, "modulo", period)
        }
        ylab <- switch(baseline, none = y,
            difference = paste("Change in", y),
            fold = paste("Fold change in", y),
            log2fold = paste("log2 fold change in", y))
        if (rolling != "none") {
            ylab <- paste0("Rolling ", rolling, " of ", ylab, " (", window,
                " samples)")
        }
        if (layout == "ridge") {
            return(.plot_ridges(df, time.col, ylab, group.by, ...))
        }
        .plot_trend(df, time.col, ylab, layout, show.individuals, smooth,
            band,
            fit, group.by, facet.by, reference, limits, xlab, ...)
    }
)

#' @importFrom utils modifyList
#' @importFrom ggplot2 ggplot aes .data geom_line geom_point geom_boxplot
#'   geom_violin geom_smooth stat_summary facet_wrap vars labs
#'   scale_x_discrete annotate geom_ribbon
#' @importFrom stats median quantile
.plot_trend <- function(df, time.col, ylab, layout, show.individuals, smooth,
        band = NULL,
        fit = NULL, group.by = NULL, facet.by = NULL, reference = NULL,
        limits = NULL, xlab = time.col, ...) {
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data$.y,
        colour = if (!is.null(group.by)) .data[[group.by]]))
    # Reference band and curve go behind the data.
    if (!is.null(limits)) {
        p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
            ymin = min(limits), ymax = max(limits), fill = "grey50",
            alpha = 0.15)
    }
    if (!is.null(reference)) {
        if (all(c("lower", "upper") %in% names(reference))) {
            p <- p + geom_ribbon(aes(x = .data[[time.col]],
                ymin = .data$lower, ymax = .data$upper), data = reference,
                inherit.aes = FALSE, fill = "grey50", alpha = 0.2)
        }
        p <- p + geom_line(aes(x = .data[[time.col]], y = .data$estimate),
            data = reference, inherit.aes = FALSE, colour = "grey30",
            linewidth = 1)
    }
    dots <- list(...)
    if (layout == "box") {
        # Individual points replace the outlier points.
        p <- p + do.call(geom_boxplot, c(list(aes(group = .data[[time.col]],
            colour = NULL)), modifyList(list(outliers = !show.individuals),
            dots)))
    }
    if (layout == "violin") {
        p <- p + do.call(geom_violin, c(list(aes(group = .data[[time.col]],
            colour = NULL)), dots))
    }
    # Individuals fade when boxes, a smooth or a fit summarise them.
    fade <- if (layout == "line" && smooth == "none" && is.null(fit) &&
            is.null(band)) list() else list(alpha = 0.3)
    if (show.individuals) {
        line_args <- if (layout == "line") modifyList(fade, dots) else fade
        p <- p + do.call(geom_line, c(list(aes(group = .data$.line)),
            line_args)) + do.call(geom_point, fade)
    }
    summary_group <- aes(group = if (is.null(group.by)) 1L else
        .data[[group.by]])
    if (!is.null(band)) {
        p <- p + stat_summary(aes(group = if (is.null(group.by)) 1L else
                .data[[group.by]], fill = if (!is.null(group.by))
                .data[[group.by]]),
            fun.min = function(v) quantile(v, band[1L], names = FALSE),
            fun.max = function(v) quantile(v, band[2L], names = FALSE),
            fun = median, geom = "ribbon", alpha = 0.2, colour = NA) +
            labs(fill = group.by)
    }
    if (smooth %in% c("mean", "median")) {
        p <- p + stat_summary(summary_group,
            fun = if (smooth == "mean") mean else median, geom = "line",
            linewidth = 1)
    } else if (smooth == "loess") {
        p <- p + geom_smooth(summary_group, method = "loess",
            formula = y ~ x, se = FALSE, linewidth = 1)
    }
    if (is.character(fit)) {
        p <- p + geom_line(aes(y = .data$.fit, group = .data$.line),
            linetype = "dashed", linewidth = 1)
    } else if (!is.null(fit)) {
        p <- p + .fit_layers(fit, time.col, group.by)
    }
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p + labs(x = xlab, y = ylab, colour = group.by)
}

# Line and interval of a results table, one per group.by group if the
# table has that column.
#' @importFrom ggplot2 aes .data geom_line geom_ribbon labs
.fit_layers <- function(fit, time.col, group.by) {
    fit_by <- if (!is.null(group.by) && group.by %in% names(fit)) group.by
    layers <- list(geom_line(aes(x = .data[[time.col]], y = .data$estimate,
            group = if (is.null(fit_by)) 1L else .data[[fit_by]],
            colour = if (!is.null(fit_by)) .data[[fit_by]]),
        data = fit, inherit.aes = FALSE, linewidth = 1))
    if (all(c("lower", "upper") %in% names(fit))) {
        layers <- c(list(geom_ribbon(aes(x = .data[[time.col]],
                ymin = .data$lower, ymax = .data$upper,
                group = if (is.null(fit_by)) 1L else .data[[fit_by]],
                fill = if (!is.null(fit_by)) .data[[fit_by]]),
            data = fit, inherit.aes = FALSE, alpha = 0.2)), layers,
            list(labs(fill = fit_by)))
    }
    layers
}

# Statistic over the last window values; NA until the window is full.
# t gives the times, used by "slope".
#' @importFrom stats var cor cov
.rolling <- function(v, stat, window, t = seq_along(v)) {
    n <- length(v)
    out <- rep(NA_real_, n)
    for (i in seq_len(n)[seq_len(n) >= window]) {
        j <- (i - window + 1L):i
        w <- v[j]
        out[i] <- switch(stat, mean = mean(w, na.rm = TRUE),
            variance = var(w, na.rm = TRUE),
            acf1 = suppressWarnings(cor(w[-1L], w[-window],
                use = "complete.obs")),
            slope = cov(t[j], w, use = "complete.obs") /
                var(t[j], na.rm = TRUE))
    }
    out
}

# Density ridges of y at each time point.
#' @importFrom ggplot2 ggplot aes .data labs
.plot_ridges <- function(df, time.col, ylab, group.by, ...) {
    df$.time <- factor(df[[time.col]])
    ggplot(df, aes(x = .data$.y, y = .data$.time,
            fill = if (!is.null(group.by)) .data[[group.by]])) +
        ggridges::geom_density_ridges(alpha = 0.6, ...) +
        labs(x = ylab, y = time.col, fill = group.by)
}
