#' Plot clusters of feature trajectories
#'
#' \code{plotProfiles()} draws one panel per cluster of features. Each
#' feature is a faint line of its mean value over the samples at each time,
#' and the mean of these lines in the cluster is a bold line.
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
#'   \code{colData(x)} that identifies subjects; it is checked for
#'   duplicated subject-time pairs. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to plot.
#'   (Default: \code{"counts"})
#'
#' @param clust.col \code{Character scalar}. Specifies the column from
#'   \code{rowData(x)} with the cluster of each feature, as written by
#'   \code{mia::addCluster()}. Features without a cluster (\code{NA}) are
#'   left out. (Default: \code{"cluster"})
#'
#' @param features \code{Character vector}. Features to plot. If
#'   \code{NULL}, plots the \code{top} features with the largest mean.
#'   (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features to plot when
#'   \code{features} is \code{NULL}. If \code{NULL}, plots all features.
#'   (Default: \code{NULL})
#'
#' @param scale \code{Logical scalar}. Should the mean trajectory of each
#'   feature be scaled to mean 0 and standard deviation 1, so that features
#'   of different abundance are compared by their shape?
#'   (Default: \code{TRUE})
#'
#' @param max.gap \code{Numeric scalar}. Breaks the lines where consecutive
#'   times are further apart than this, in \code{time.col} units. For
#'   visits, \code{NULL} breaks lines at visits without a sample; otherwise
#'   \code{NULL} never breaks lines. (Default: \code{NULL})
#'
#' @param events \code{data.frame}. Events or intervals to mark in every
#'   panel, one per row, with columns \code{start} and optionally
#'   \code{end} and \code{label}, in the units of \code{time.col} (visit
#'   levels for visits). Rows with an \code{end} are drawn as grey bands
#'   with their label at the top, the others as vertical lines.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}} for the feature lines.
#'
#' @details
#' The value of a feature at a time is its mean over all samples at that
#' time, over subjects and groups alike; with irregular sampling, where
#' most times have a single sample, give the times as a factor of visits.
#' The clusters come from upstream, such as a clustering of the features
#' by \code{mia::addCluster()} on a scaled assay; they are drawn in the
#' order of the levels of \code{clust.col}, or sorted. The panel titles
#' give each cluster and its number of plotted features.
#'
#' The panels follow the trajectory clusters of Bodein et al. (2019,
#' doi:10.3389/fgene.2019.00963, Figs 6 and 7), where each line is the
#' scaled profile of a feature, and the clusters of taxa responding to a
#' diet in David et al. (2014, doi:10.1186/gb-2014-15-7-r89, Fig 3).
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Bodein A, et al. (2019). A generic multivariate framework for the
#' integration of microbiome longitudinal studies with other data types.
#' \emph{Frontiers in Genetics} 10:963.
#'
#' David LA, et al. (2014). Host lifestyle affects human microbiota on
#' daily timescales. \emph{Genome Biology} 15(7):R89.
#'
#' @name plotProfiles
#'
#' @examples
#' # 16 species in three bioreactors, clustered by the shape of their
#' # profiles: hierarchical clustering of the abundances scaled per
#' # species, cut into four clusters
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "relabundance")
#' z <- t(scale(t(SummarizedExperiment::assay(tse, "relabundance"))))
#' SummarizedExperiment::rowData(tse)$cluster <- factor(
#'     cutree(hclust(dist(z), method = "ward.D2"), k = 4))
#' plotProfiles(tse, time.col = "Time_hr_num", group = "StudyIdentifier",
#'     assay.type = "relabundance")
#'
#' # Unscaled relative abundances, with the batch phase and the addition
#' # of Blautia hydrogenotrophica at 152 h marked
#' events <- data.frame(start = c(0, 152), end = c(24, NA),
#'     label = c("batch", "Blautia added"))
#' plotProfiles(tse, time.col = "Time_hr_num", group = "StudyIdentifier",
#'     assay.type = "relabundance", scale = FALSE, events = events)
NULL

#' @rdname plotProfiles
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom stats aggregate ave na.pass
#' @importFrom methods setMethod
#' @export
setMethod("plotProfiles", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            clust.col = "cluster", features = NULL, top = NULL,
            scale = TRUE, max.gap = NULL, events = NULL, ...) {
        .require_column(x, clust.col, where = "rowData")
        .check_flag(scale)
        x <- x[!is.na(rowData(x)[[clust.col]]), ]
        df <- .feature_table(x, time.col, group, assay.type, features,
            top %||% Inf, level = "group")
        events <- .check_events(events, x, time.col)
        time <- colData(x)[[time.col]]

        cluster <- rowData(x)[[clust.col]][match(df$.feature, rownames(x))]
        lv <- .levels_of(cluster)
        df$.cluster <- factor(as.character(cluster), lv)
        n <- lengths(lapply(split(df$.feature, df$.cluster), unique))
        levels(df$.cluster) <- paste0(lv, " (n = ", n, ")")
        if (scale) {
            df$.value <- ave(df$.value, df$.feature, FUN = .zscore)
        }
        means <- aggregate(.value ~ .cluster + .time, data = df,
            FUN = mean, na.action = na.pass)

        # Both tables have the same times, so the same gaps.
        visit <- is.factor(time)
        df <- .break_gaps(.order_time(df, ".time", ".feature"), ".time",
            ".feature", max.gap, visit)
        means <- .break_gaps(.order_time(means, ".time", ".cluster"),
            ".time", ".cluster", max.gap, visit)
        df$.time <- .time_restore(df$.time, time)
        means$.time <- .time_restore(means$.time, time)
        ylab <- paste(if (scale) "Scaled mean" else "Mean", assay.type)
        .plot_profiles(df, means, time.col, ylab, events, levels(time), ...)
    }
)

# v minus its mean, divided by its standard deviation; a constant v
# becomes 0.
#' @importFrom stats sd
.zscore <- function(v) {
    s <- sd(v, na.rm = TRUE)
    (v - mean(v, na.rm = TRUE)) / if (isTRUE(s > 0)) s else 1
}

#' @importFrom utils modifyList
#' @importFrom ggplot2 ggplot aes .data geom_line facet_wrap vars labs
#'   scale_x_continuous
.plot_profiles <- function(df, means, time.col, ylab, events, levels,
        ...) {
    ggplot(df, aes(x = .data$.time, y = .data$.value)) +
        .events_layer(events, "bands", text = TRUE, x = df$.time) +
        do.call(geom_line, c(list(aes(group = interaction(.data$.feature,
            .data$.segment))), modifyList(list(colour = "grey50",
            alpha = 0.5), list(...)))) +
        geom_line(aes(group = .data$.segment), data = means,
            colour = .okabe_ito[[1L]], linewidth = 1.2) +
        list(.events_layer(events, "markers"),
            facet_wrap(vars(.data$.cluster)),
            if (!is.null(levels)) scale_x_continuous(
                breaks = seq_along(levels), labels = levels),
            labs(x = time.col, y = ylab),
            .theme_tempoPlot())
}
