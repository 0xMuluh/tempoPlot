#' Plot a lag profile: a variogram or autocorrelations
#'
#' \code{plotLagProfile()} draws a statistic against time lag: the
#' variogram \eqn{D(L)}, how far samples of the same subject are apart at
#' each lag, with its confidence interval, as stored by
#' \code{miaTime::addVariogram()}; or the autocorrelation of each feature
#' in each subject.
#'
#' @param x a \code{data.frame} with a variogram or an autocorrelation
#'   table, or a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object with such a table in its \code{metadata}.
#'
#' @param name \code{Character scalar}. Specifies the entry of
#'   \code{metadata(x)} holding the table. (Default: \code{"variogram"})
#'
#' @param type \code{Character scalar}. For variograms, draws the
#'   dissimilarity at each lag (\code{"dissimilarity"}) or, for a variogram
#'   computed with \code{split.by}, the ratio between its levels
#'   (\code{"ratio"}), on a log scale with a line at 1.
#'   (Default: \code{"dissimilarity"})
#'
#' @param feature.col \code{Character scalar}. For autocorrelations,
#'   specifies the column of feature names; \code{NULL} for a table of one
#'   series. (Default: \code{"FeatureID"})
#'
#' @param group \code{Character scalar}. For autocorrelations, specifies
#'   the column that identifies subjects. (Default: \code{NULL})
#'
#' @param colour.by \code{Character scalar}. For autocorrelations,
#'   specifies the column that colours the lines, such as
#'   \code{feature.col}. If \code{NULL}, all lines take one colour.
#'   (Default: \code{NULL})
#'
#' @param facet.by \code{Character scalar}. Specifies the column that
#'   splits the table into panels. (Default: \code{group})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_linerange]{geom_pointrange}} for variograms
#'   and \code{\link[ggplot2:geom_path]{geom_line}} for autocorrelations.
#'
#' @details
#' The variogram is a table with the columns \code{lag}, \code{estimate},
#' \code{lower} and \code{upper}, and optionally \code{type} and, between
#' \code{lag} and \code{type}, the \code{split.by} column. With
#' \code{split.by}, each level or comparison is drawn in its own colour,
#' slightly shifted apart at each lag.
#'
#' An autocorrelation table has the columns \code{lag} and \code{acf},
#' one row per feature, subject and lag, as in Gibbons et al. (2017,
#' doi:10.1371/journal.pcbi.1005364, Fig. 5). Each feature in each subject
#' is a line, by default in one panel per subject. A column
#' \code{time.lag} (the lag in time units) replaces \code{lag} on the x
#' axis. With a column \code{n.pairs}, the pairs of samples at each lag,
#' a band shows \eqn{\pm 1.96 / \sqrt{n}}{+-1.96 / sqrt(n)}, the range of
#' autocorrelations expected without dependence between samples.
#'
#' Lags without an estimate are dropped.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotLagProfile
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- mia::transformAssay(hitchip1006, method = "hellinger")
#'
#' # Squared Hellinger distance between samples of a subject L time units
#' # apart
#' tse <- miaTime::addVariogram(tse, time.col = "time", group = "subject",
#'     assay.type = "hellinger", method = "euclidean", squared = TRUE,
#'     time.lag = c(1.5, 2, 2.5, 3, 4, 6, 8, 9))
#' plotLagProfile(tse)
#'
#' # The same for women and men
#' tse <- miaTime::addVariogram(tse, time.col = "time", group = "subject",
#'     assay.type = "hellinger", method = "euclidean", squared = TRUE,
#'     time.lag = c(2, 6, 8), split.by = "sex", name = "variogram_sex")
#' plotLagProfile(tse, name = "variogram_sex")
#'
#' # Their ratio
#' plotLagProfile(tse, name = "variogram_sex", type = "ratio")
#'
#' # Hourly samples with reads of four vessels of an artificial gut, hours
#' # 476 to 592, one sample per hour, and the five most abundant taxa
#' data("SilvermanAGutData", package = "miaTime")
#' tse <- SilvermanAGutData
#' tse <- tse[, tse$SampleType == "Hourly" &
#'     tse$Normal_Noise_Sample == "Normal" & tse$Time.hr %in% 476:592 &
#'     colSums(SummarizedExperiment::assay(tse)) > 0]
#' tse <- tse[, !duplicated(paste(tse$Vessel, tse$Time.hr))]
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' tse <- tse[mia::getTop(tse, top = 5), ]
#'
#' # Illustrative autocorrelation table in the form of miaTime's
#' # getAutocorrelation(), from stats::acf() on the hourly grid; hours
#' # without a sample are NA
#' acfs <- do.call(rbind, lapply(unique(tse$Vessel), function(v) {
#'     j <- match(476:592, tse$Time.hr[tse$Vessel == v])
#'     m <- SummarizedExperiment::assay(tse, "relabundance")[,
#'         tse$Vessel == v][, j]
#'     present <- !is.na(j)
#'     do.call(rbind, lapply(rownames(m), function(feature) {
#'         r <- stats::acf(m[feature, ], lag.max = 12,
#'             na.action = stats::na.pass, plot = FALSE)
#'         data.frame(FeatureID = feature, Vessel = v, lag = 1:12,
#'             acf = r$acf[-1L], n.pairs = vapply(1:12, function(l)
#'                 sum(present[-seq_len(l)] & head(present, -l)),
#'                 integer(1)))
#'     }))
#' }))
#' plotLagProfile(acfs, group = "Vessel", colour.by = "FeatureID")
NULL

#' @rdname plotLagProfile
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotLagProfile", signature = c(x = "SummarizedExperiment"),
    function(x, name = "variogram", ...) {
        if (!.is_one_name(name) || is.null(metadata(x)[[name]])) {
            stop("'name' must be an entry of metadata(x); compute it with ",
                "miaTime::addVariogram().", call. = FALSE)
        }
        plotLagProfile(as.data.frame(metadata(x)[[name]]), ...)
    }
)

#' @rdname plotLagProfile
#' @importFrom methods setMethod
#' @importFrom ggplot2 facet_wrap vars
#' @export
setMethod("plotLagProfile", signature = c(x = "data.frame"),
    function(x, type = c("dissimilarity", "ratio"),
            feature.col = "FeatureID", group = NULL, colour.by = NULL,
            facet.by = group, ...) {
        type <- match.arg(type)
        # Only facet.by applies to variograms.
        acf <- "acf" %in% names(x)
        .check_table_columns(x, feature.col = if (acf) feature.col,
            group = if (acf) group, colour.by = if (acf) colour.by,
            facet.by = facet.by)
        p <- if (acf) {
            .plot_autocorrelation(x, feature.col, group, colour.by,
                facet.by, ...)
        } else {
            .variogram_profile(x, type, ...)
        }
        if (!is.null(facet.by)) {
            p <- p + facet_wrap(vars(.data[[facet.by]]))
        }
        p
    }
)

# Plot of the variogram table x: estimates of type with their intervals.
.variogram_profile <- function(x, type, ...) {
    if (!all(c("lag", "estimate", "lower", "upper") %in% names(x))) {
        stop("'x' must have the columns lag and acf (autocorrelations), ",
            "or lag, estimate, lower and upper (a variogram).",
            call. = FALSE)
    }
    if ("type" %in% names(x)) {
        x <- x[x$type == type, , drop = FALSE]
    } else if (type == "ratio") {
        x <- x[0L, , drop = FALSE]
    }
    x <- x[!is.na(x$estimate), , drop = FALSE]
    if (nrow(x) == 0L) {
        stop("'x' has no ", type, " estimates",
            if (type == "ratio") "; compute the variogram with 'split.by'",
            ".", call. = FALSE)
    }
    # The split.by column sits between lag and type; its ratio rows
    # name the comparisons, such as "male / female".
    split <- setdiff(names(x)[seq_len(match("type", names(x), 0L))],
        c("lag", "type"))
    .plot_lag_profile(x, if (length(split) == 1L) split, type, ...)
}

#' @importFrom ggplot2 ggplot aes .data geom_line geom_pointrange labs
#'   position_dodge geom_hline scale_y_log10
.plot_lag_profile <- function(df, split, type = "dissimilarity", ...) {
    # Levels at the same lag are shifted apart so their intervals show.
    dodge <- position_dodge(width = 0.02 * diff(range(df$lag)))
    p <- ggplot(df, aes(x = .data$lag, y = .data$estimate,
            colour = if (!is.null(split)) .data[[split]],
            group = if (!is.null(split)) .data[[split]] else 1L)) +
        geom_line(position = dodge) +
        geom_pointrange(aes(ymin = .data$lower, ymax = .data$upper),
            position = dodge, ...) +
        .scale_categorical("colour") +
        labs(x = "lag", y = if (type == "ratio") "Ratio of D(L)" else "D(L)",
            colour = split)
    if (type == "ratio") {
        # A ratio of 1 means the same dissimilarity in both levels.
        p <- p + geom_hline(yintercept = 1, linetype = "dashed",
                colour = "grey40") +
            scale_y_log10()
    }
    p + .theme_tempoPlot()
}

# Lines of the autocorrelation table x against lag (time.lag if given),
# one per feature and subject, with a band of +-1.96 / sqrt(n.pairs)
# around 0 for each subject if given.
#' @importFrom ggplot2 ggplot aes .data geom_hline geom_line geom_ribbon
#'   labs
#' @importFrom utils modifyList
.plot_autocorrelation <- function(x, feature.col, group, colour.by,
        facet.by, ...) {
    x$.lag <- if ("time.lag" %in% names(x)) x$time.lag else x$lag
    x <- x[!is.na(x$acf), , drop = FALSE]
    ids <- c(feature.col, group)
    x$.line <- if (length(ids)) interaction(x[ids], drop = TRUE) else 1L
    # Without colour.by, all lines take the first colour.
    fixed <- if (is.null(colour.by)) list(colour = .okabe_ito[[1L]])
    p <- ggplot(x, aes(x = .data$.lag, y = .data$acf)) +
        geom_hline(yintercept = 0, colour = "grey40")
    if ("n.pairs" %in% names(x)) {
        band <- unique(x[x$n.pairs > 0,
            unique(c(group, facet.by, ".lag", "n.pairs")), drop = FALSE])
        band$.band <- 1.96 / sqrt(band$n.pairs)
        p <- p + geom_ribbon(aes(x = .data$.lag, ymin = -.data$.band,
                ymax = .data$.band,
                group = if (!is.null(group)) .data[[group]]), data = band,
                inherit.aes = FALSE, fill = "grey50", alpha = 0.2)
    }
    p + do.call(geom_line, c(list(aes(group = .data$.line,
            colour = if (!is.null(colour.by)) .data[[colour.by]])),
            modifyList(as.list(fixed), list(...)))) +
        .scale_by(x, colour.by, "colour") +
        labs(x = "lag", y = "autocorrelation", colour = colour.by) +
        .theme_tempoPlot()
}
