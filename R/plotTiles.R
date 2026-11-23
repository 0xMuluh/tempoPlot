#' Plot features over time as tiles
#'
#' \code{plotTiles()} draws a features by time heatmap of an assay, one
#' panel per subject or per group of subjects.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object, or a \code{data.frame} of results with one row per feature
#'   and time point, such as differential abundance at each visit.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)}, or of a \code{data.frame} \code{x}, that contains
#'   times: numeric, \code{Date}, \code{POSIXct} or a factor of visits.
#'   Time differences, as in \code{max.gap}, are in days for \code{Date}
#'   and hours for \code{POSIXct} times.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)}, or of a \code{data.frame} \code{x}, that identifies
#'   subjects. If \code{NULL}, all samples, or all rows, form one series.
#'   (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to plot.
#'   (Default: \code{"counts"})
#'
#' @param features \code{Character vector}. Features to plot. If
#'   \code{NULL}, plots the \code{top} features with the largest mean, or
#'   the largest mean absolute effect for a \code{data.frame}.
#'   (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features to plot when
#'   \code{features} is \code{NULL}. (Default: \code{20})
#'
#' @param value \code{Character scalar}. Plots the assay values
#'   (\code{"raw"}), their change from the first sample of each subject
#'   (\code{"change"}), or presence, a value above \code{detection}
#'   (\code{"occupancy"}). Effects from a results table are drawn by the
#'   \code{data.frame} method. (Default: \code{"raw"})
#'
#' @param detection \code{Numeric scalar}. With \code{value =
#'   "occupancy"}, a feature is present in a sample when its value is above
#'   this. (Default: \code{0})
#'
#' @param effect \code{Character scalar}. For a \code{data.frame}, the
#'   column of effects, such as log fold changes at each time point.
#'
#' @param p.col \code{Character scalar}. For a \code{data.frame}, the column
#'   of p-values or adjusted p-values; cells below \code{threshold} are
#'   marked with a dot. (Default: \code{NULL})
#'
#' @param threshold \code{Numeric scalar}. Significance level for
#'   \code{p.col}. (Default: \code{0.05})
#'
#' @param feature.col \code{Character scalar}. For a \code{data.frame}, the
#'   column of feature names. (Default: \code{"feature"})
#'
#' @param order.by \code{Character scalar}. Orders features by their mean
#'   plotted value, largest at the top (\code{"mean"}), or by hierarchical
#'   clustering of their profiles (\code{"cluster"}). (Default: \code{"mean"})
#'
#' @param level \code{Character scalar}. Draws one panel per subject
#'   (\code{"subject"}) or averages the subjects of each \code{group.by}
#'   group at each time (\code{"group"}); occupancy then becomes the share
#'   of subjects with the feature. (Default: \code{"subject"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that groups subjects when \code{level = "group"}. If
#'   \code{NULL}, all subjects form one group. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Leaves blank the gaps between
#'   consecutive times that are further apart than this, in \code{time.col}
#'   units: a tile reaches at most \code{max.gap / 2} from its time.
#'   \code{NULL} never leaves gaps. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_tile]{geom_tile}}.
#'
#' @details
#' Each sample covers the time from halfway to its previous sample to
#' halfway to its next one, so the time axis keeps the spacing of the
#' sampling times. Visits are evenly spaced. Missing values are drawn
#' grey, times without a sample blank.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotTiles
#'
#' @examples
#' # 16 species in three bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "relabundance")
#' plotTiles(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     assay.type = "relabundance")
#'
#' # Change from the start, clustered, averaged over the bioreactors
#' plotTiles(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     assay.type = "relabundance", value = "change",
#'     order.by = "cluster", level = "group")
#'
#' # Clr difference between men and women at each visit, with Wilcoxon
#' # q-values below 0.05 marked
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006[, !is.na(hitchip1006$sex)]
#' tse$visit <- ave(tse$time, tse$subject, FUN = rank)
#' tse <- mia::transformAssay(tse, method = "clr", pseudocount = 1)
#' mat <- SummarizedExperiment::assay(tse, "clr")
#' res <- do.call(rbind, lapply(1:3, function(v) {
#'     m <- mat[, tse$visit == v]
#'     male <- tse$sex[tse$visit == v] == "male"
#'     data.frame(feature = rownames(m), visit = v,
#'         difference = rowMeans(m[, male]) - rowMeans(m[, !male]),
#'         p = apply(m, 1, function(f) wilcox.test(f[male],
#'             f[!male])$p.value))
#' }))
#' res$q <- p.adjust(res$p, method = "fdr")
#' plotTiles(res, time.col = "visit", effect = "difference", p.col = "q",
#'     top = 15)
#'
#' # Rate of change of each species in each bioreactor, from miaTime
#' tse <- mia::transformAssay(minimalgut, method = "relabundance")
#' change <- miaTime::getShortTermChange(tse, assay.type = "relabundance",
#'     time.col = "Time.hr", group = "StudyIdentifier")
#' plotTiles(as.data.frame(change), time.col = "Time.hr",
#'     effect = "rate_of_change", feature.col = "FeatureID",
#'     group = "StudyIdentifier", top = 8)
NULL

#' @rdname plotTiles
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotTiles", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            features = NULL, top = 20,
            value = c("raw", "change", "occupancy"), detection = 0,
            order.by = c("mean", "cluster"),
            level = c("subject", "group"), group.by = NULL, max.gap = NULL,
            ...) {
        value <- match.arg(value)
        order.by <- match.arg(order.by)
        level <- match.arg(level)
        df <- .feature_table(x, time.col, group, assay.type, features, top,
            value, level, group.by, detection)
        time <- colData(x)[[time.col]]
        df$.feature <- factor(df$.feature,
            levels = .order_features(df, order.by))
        df <- .add_tile_edges(df, time, max.gap)
        panel.col <- if (level == "subject") group else group.by
        .plot_tiles(df, time.col, assay.type, value, level, panel.col,
            levels(time), ...)
    }
)

#' @rdname plotTiles
#' @importFrom methods setMethod
#' @export
setMethod("plotTiles", signature = c(x = "data.frame"),
    function(x, time.col, effect, p.col = NULL, threshold = 0.05,
            features = NULL, top = 20, order.by = c("mean", "cluster"),
            feature.col = "feature", group = NULL, max.gap = NULL, ...) {
        order.by <- match.arg(order.by)
        .check_effects(x, time.col, effect, p.col, threshold, feature.col,
            group)
        time <- x[[time.col]]
        if (is.factor(time)) {
            time <- droplevels(time)
        }
        # Features with the largest mean absolute effect.
        size <- tapply(abs(x[[effect]]), x[[feature.col]], mean,
            na.rm = TRUE)
        mat <- matrix(size, dimnames = list(names(size), NULL))
        features <- .select_features(mat, features, top)
        df <- data.frame(.feature = as.character(x[[feature.col]]),
            .panel = if (is.null(group)) "" else as.character(x[[group]]),
            .time = .time_numeric(time),
            .value = x[[effect]])
        if (!is.null(p.col)) {
            df$.significant <- !is.na(x[[p.col]]) & x[[p.col]] < threshold
        }
        df <- df[df$.feature %in% features, , drop = FALSE]
        if (anyDuplicated(df[c(".panel", ".feature", ".time")])) {
            stop("'x' must have one row per feature and time point",
                if (!is.null(group)) " of each subject", ".", call. = FALSE)
        }
        df$.feature <- factor(df$.feature,
            levels = .order_features(df, order.by))
        df <- .add_tile_edges(df, time, max.gap)
        p <- .plot_tiles(df, time.col, effect, "effect", "subject", group,
            levels(time), ...)
        if (!is.null(p.col)) {
            p <- p + geom_point(aes(x = .data$.mid, y = .data$.feature),
                data = df[df$.significant, , drop = FALSE], size = 1,
                inherit.aes = FALSE) +
                labs(caption = paste("Dots:", p.col, "<", threshold))
        }
        p
    }
)

# Check the columns, times and threshold of a results table x for
# plotTiles().
.check_effects <- function(x, time.col, effect, p.col, threshold,
        feature.col, group = NULL) {
    for (arg in c("time.col", "effect", "p.col", "feature.col", "group")) {
        v <- get(arg)
        if (arg %in% c("p.col", "group") && is.null(v)) next
        if (!.is_one_name(v) || !v %in% names(x)) {
            stop("'", arg, "' must be one column name of 'x'.",
                call. = FALSE)
        }
    }
    if (!is.numeric(x[[effect]])) {
        stop("'effect' must be a numeric column.", call. = FALSE)
    }
    .check_number(threshold, 0, 1, open = TRUE)
    .time_type(x[[time.col]])
}

# Feature levels, bottom to top: by mean value (largest on top) or by the
# order of a hierarchical clustering of the feature profiles.
#' @importFrom stats dist hclust
.order_features <- function(df, order.by) {
    if (order.by == "mean") {
        m <- tapply(df$.value, df$.feature, mean, na.rm = TRUE)
        return(names(sort(m, na.last = FALSE)))
    }
    key <- paste(df$.panel, df$.time)
    wide <- tapply(df$.value, list(df$.feature, key), mean)
    wide[is.na(wide)] <- 0
    if (nrow(wide) < 2L) {
        return(rownames(wide))
    }
    rownames(wide)[hclust(dist(wide))$order]
}

# Tile centre and width of each time within each panel, for the times of
# time.col. A tile reaches at most max.gap / 2 from its time, so longer
# gaps stay blank. Centres are returned in the class of time.
.add_tile_edges <- function(df, time, max.gap = NULL) {
    .check_number(max.gap, 0, open = TRUE, null = TRUE)
    half <- (max.gap %||% Inf) / 2
    edges <- lapply(split(df$.time, df$.panel), function(t) {
        t <- sort(unique(t))
        e <- .tile_edges(t, is.factor(time))
        data.frame(.time = t, .lower = pmax(e$lower, t - half),
            .upper = pmin(e$upper, t + half))
    })
    edges <- do.call(rbind, Map(cbind, .panel = names(edges), edges))
    df <- merge(df, edges, by = c(".panel", ".time"), sort = FALSE)
    df$.mid <- (df$.lower + df$.upper) / 2
    df$.width <- df$.upper - df$.lower
    df$.mid <- .time_restore(df$.mid, time)
    df
}

#' @importFrom ggplot2 ggplot aes .data geom_tile geom_point labs facet_wrap
#'   vars scale_x_continuous
.plot_tiles <- function(df, time.col, assay.type, value, level, panel.col,
        levels, ...) {
    lab <- switch(value, raw = assay.type,
        change = paste(assay.type, "change"), effect = assay.type,
        occupancy = if (level == "group") "occupancy" else "present")
    p <- ggplot(df, aes(x = .data$.mid, y = .data$.feature)) +
        geom_tile(aes(width = .data$.width, fill = .data$.value), ...) +
        switch(value, raw = .scale_sequential("fill"),
            occupancy = .scale_sequential("fill", limits = c(0, 1)),
            .scale_diverging("fill")) +
        labs(x = time.col, y = NULL, fill = lab) +
        .theme_tempoPlot()
    if (!is.null(levels)) {
        p <- p + scale_x_continuous(breaks = seq_along(levels),
            labels = levels)
    }
    if (!is.null(panel.col)) {
        p <- p + facet_wrap(vars(.data$.panel))
    }
    p
}
