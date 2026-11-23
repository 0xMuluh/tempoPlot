#' Plot the timeline of one subject
#'
#' \code{plotTimeline()} stacks tracks of one subject on a shared time
#' axis: the composition of the most abundant features, numeric variables
#' as lines, a strip of states and a strip of events.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits. Time differences, as in
#'   \code{max.gap}, are in days for \code{Date} and hours for
#'   \code{POSIXct} times.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   belong to one subject. (Default: \code{NULL})
#'
#' @param subject \code{Character scalar}. The subject to draw, a value of
#'   the \code{group} column. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of
#'   abundances, also read for \code{y} that are features.
#'   (Default: \code{"counts"})
#'
#' @param y \code{Character vector}. Specifies the variables drawn as
#'   lines, one track each, as in \code{\link{plotTrend}}: columns of
#'   \code{colData(x)}, rows of \code{assay(x, assay.type)} or
#'   \code{"<dimred>:<component>"}. (Default: \code{NULL})
#'
#' @param state.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains the state of each sample, such as a
#'   community type, drawn as a strip of tiles. (Default: \code{NULL})
#'
#' @param events \code{data.frame}. Events to draw in a strip at the
#'   bottom, one per row, with columns \code{start}, and optionally
#'   \code{end}, \code{label} and the \code{group} column. \code{start} and
#'   \code{end} have the type of \code{time.col}; for visits, they are
#'   visit levels. Events with an \code{end} are drawn as bands, the others
#'   as vertical lines; \code{label} sets their fill or line type. Events
#'   with a \code{group} column are drawn for their subject only.
#'   (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Breaks lines, and areas of the
#'   composition, where consecutive samples are further apart than this,
#'   in \code{time.col} units. For visits, \code{NULL} breaks them at
#'   visits without a sample; otherwise \code{NULL} never breaks them.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to \code{\link{plotComposition}},
#'   such as \code{features}, \code{top}, \code{as.relative} or
#'   \code{layout}.
#'
#' @details
#' The tracks share the x axis: the same limits, from half the spacing of
#' samples before the first to half after the last, widened to include
#' the events, and the same breaks. Only the bottom track shows the time
#' axis, and the legends of all tracks are collected on the right. The
#' layout follows the subject timelines of Gajer et al. (2012), Figs 3 and
#' 5, Zhou et al. (2019), Fig. 6b, and David et al. (2014), Fig. 1.
#'
#' Draws with the \pkg{patchwork} package. A single track is changed with
#' \code{p[[i]] <- p[[i]] + ...}, all of them with \code{p & ...}.
#'
#' @return
#' A \code{patchwork} object, not a single \code{ggplot}: the tracks from
#' top to bottom.
#'
#' @references
#' David LA, et al. (2014). Host lifestyle affects human microbiota on
#' daily timescales. \emph{Genome Biology} 15(7):R89.
#' doi:10.1186/gb-2014-15-7-r89
#'
#' Gajer P, et al. (2012). Temporal dynamics of the human vaginal
#' microbiota. \emph{Science Translational Medicine} 4(132):132ra52.
#' doi:10.1126/scitranslmed.3003605
#'
#' Zhou W, et al. (2019). Longitudinal multi-omics of host-microbe
#' dynamics in prediabetes. \emph{Nature} 569(7758):663-671.
#' doi:10.1038/s41586-019-1236-x
#'
#' @name plotTimeline
#'
#' @examples
#' # A bioreactor of a defined community: species, butyrate, the dominant
#' # species of each sample and the perturbations of its schedule
#' data("minimalgut", package = "miaTime")
#' tse <- mia::addDominant(minimalgut, name = "dominant")
#' events <- data.frame(
#'     start = c(0, 152, 432),
#'     end = c(24, NA, 436),
#'     label = c("batch", "Blautia added", "double dilution rate"))
#' plotTimeline(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     subject = "Bioreactor A", y = "Butyrate", state.col = "dominant",
#'     events = events, top = 8, layout = "area")
NULL

#' @rdname plotTimeline
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 theme element_blank
#' @importFrom methods setMethod
#' @export
setMethod("plotTimeline", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, subject = NULL,
            assay.type = "counts", y = NULL, state.col = NULL,
            events = NULL, max.gap = NULL, ...) {
        if (!requireNamespace("patchwork", quietly = TRUE)) {
            stop("plotTimeline() requires the patchwork package.",
                call. = FALSE)
        }
        .check_column(x, group)
        if (!is.null(group)) {
            if (!(length(subject) == 1L &&
                    subject %in% colData(x)[[group]])) {
                stop("'subject' must be one value of the 'group' column.",
                    call. = FALSE)
            }
            x <- x[, colData(x)[[group]] %in% subject]
        }
        .check_longitudinal(x, time.col)
        .check_column(x, state.col)
        events <- .check_events(events, x, time.col, group)
        if (!is.null(group) && group %in% names(events)) {
            events <- events[events[[group]] %in% subject, , drop = FALSE]
        }
        time <- colData(x)[[time.col]]

        tracks <- c(
            list(plotComposition(x, time.col, assay.type = assay.type,
                max.gap = max.gap, ...)),
            lapply(y, function(v) plotTrend(x, v, time.col,
                assay.type = assay.type, max.gap = max.gap)),
            if (!is.null(state.col)) list(.plot_state_strip(x, time.col,
                state.col) + .scale_time(time)),
            if (NROW(events)) list(.plot_event_strip(events, time.col) +
                .scale_time(time)))
        xlim <- .timeline_limits(time, events)
        n <- length(tracks)
        tracks <- lapply(seq_len(n), function(i) tracks[[i]] + xlim +
            if (i < n) theme(axis.title.x = element_blank(),
                axis.text.x = element_blank()))
        heights <- c(3, rep(1.5, length(y)), rep(0.4, n - 1L - length(y)))
        patchwork::wrap_plots(tracks, ncol = 1L, heights = heights,
            guides = "collect") & theme(legend.position = "right")
    }
)

# x axis of every track: the sample times with half their spacing at both
# ends (all visits), widened to the events.
#' @importFrom ggplot2 coord_cartesian
.timeline_limits <- function(time, events) {
    visit <- is.factor(time)
    t <- if (visit) seq_along(levels(time)) else
        sort(unique(.time_numeric(time)))
    edges <- .tile_edges(t, visit)
    lim <- range(edges$lower, edges$upper, .time_numeric(events$start),
        .time_numeric(events$end), na.rm = TRUE)
    coord_cartesian(xlim = .time_restore(lim, time),
        expand = c(left = FALSE, right = FALSE))
}

# Discrete x scale with all visits for the strips, which may hold no
# samples at some visits (states) or no factor at all (events).
#' @importFrom ggplot2 scale_x_discrete
.scale_time <- function(time) {
    if (is.factor(time)) scale_x_discrete(limits = levels(time))
}

#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes .data geom_tile labs
.plot_state_strip <- function(x, time.col, state.col) {
    df <- as.data.frame(colData(x)[, c(time.col, state.col)],
        optional = TRUE)
    state <- df[[state.col]]
    df$.state <- factor(as.character(state), levels = .levels_of(state))
    ggplot(df, aes(x = .data[[time.col]], y = state.col,
            fill = .data$.state)) +
        geom_tile() +
        .scale_categorical("fill") +
        labs(x = time.col, y = NULL, fill = state.col) +
        .theme_tempoPlot()
}

#' @importFrom ggplot2 ggplot labs scale_y_continuous
.plot_event_strip <- function(events, time.col) {
    ggplot() +
        .events_layer(events, "bands") +
        .events_layer(events, "markers") +
        scale_y_continuous(limits = c(0, 1), breaks = 0.5,
            labels = "events") +
        .scale_categorical("fill") +
        labs(x = time.col, y = NULL) +
        .theme_tempoPlot()
}
