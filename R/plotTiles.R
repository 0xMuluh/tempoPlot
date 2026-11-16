#' Plot features over time as tiles
#'
#' \code{plotTiles()} draws a features by time heatmap of an assay, one
#' panel per subject or per group of subjects.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to plot.
#'   (Default: \code{"counts"})
#'
#' @param features \code{Character vector}. Features to plot. If
#'   \code{NULL}, plots the \code{top} features with the largest mean.
#'   (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features to plot when
#'   \code{features} is \code{NULL}. (Default: \code{20})
#'
#' @param value \code{Character scalar}. Plots the assay values
#'   (\code{"raw"}), their change from the first sample of each subject
#'   (\code{"change"}), or presence, a value above zero
#'   (\code{"occupancy"}). (Default: \code{"raw"})
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
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_tile]{geom_tile}}.
#'
#' @details
#' Each sample covers the time from halfway to its previous sample to
#' halfway to its next one, so the time axis keeps the spacing of the
#' sampling times. Visits are evenly spaced. Missing values are drawn
#' blank.
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
NULL

#' @rdname plotTiles
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom methods setMethod
#' @export
setMethod("plotTiles", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            features = NULL, top = 20,
            value = c("raw", "change", "occupancy"),
            order.by = c("mean", "cluster"),
            level = c("subject", "group"), group.by = NULL, ...) {
        value <- match.arg(value)
        order.by <- match.arg(order.by)
        level <- match.arg(level)
        .check_longitudinal(x, time.col, group)
        type <- .time_type(x, time.col)
        .check_column(x, group.by)
        if (!.is_one_name(assay.type) || !assay.type %in% assayNames(x)) {
            stop("'assay.type' must be one of assayNames(x).",
                call. = FALSE)
        }
        mat <- as.matrix(assay(x, assay.type))
        features <- .select_features(mat, features, top)
        mat <- mat[features, , drop = FALSE]

        time <- colData(x)[[time.col]]
        subject <- if (is.null(group)) rep("", ncol(x)) else
            as.character(colData(x)[[group]])
        if (value == "change") {
            first <- vapply(split(seq_len(ncol(x)), subject),
                function(i) i[which.min(.time_numeric(time[i]))],
                integer(1))
            mat <- mat - mat[, first[match(subject, names(first))],
                drop = FALSE]
        } else if (value == "occupancy") {
            mat <- (mat > 0) + 0
        }

        panel <- if (level == "subject") subject else if
            (is.null(group.by)) rep("", ncol(x)) else
            as.character(colData(x)[[group.by]])
        df <- data.frame(
            .feature = rep(features, times = ncol(x)),
            .panel = rep(panel, each = length(features)),
            .time = rep(.time_numeric(time), each = length(features)),
            .value = as.vector(mat))
        if (level == "group") {
            df <- aggregate(.value ~ .feature + .panel + .time, data = df,
                FUN = mean, na.action = na.pass)
        }
        df$.feature <- factor(df$.feature,
            levels = .order_features(df, order.by))
        df <- .add_tile_edges(df, visit = type == "visit")
        if (type == "date") {
            df$.mid <- as.Date(df$.mid, origin = "1970-01-01")
        }
        levels <- if (type == "visit") levels(time)
        panel.col <- if (level == "subject") group else group.by
        .plot_tiles(df, time.col, assay.type, value, level, panel.col,
            levels, ...)
    }
)

# Features given, checked against rownames, or the top ones by mean.
.select_features <- function(mat, features, top) {
    if (!is.null(features)) {
        if (!is.character(features) || !all(features %in% rownames(mat))) {
            stop("'features' must be row names of 'x'.", call. = FALSE)
        }
        return(unique(features))
    }
    if (!(is.numeric(top) && length(top) == 1L && top >= 1)) {
        stop("'top' must be a single number of at least 1.", call. = FALSE)
    }
    m <- rowMeans(mat, na.rm = TRUE)
    rownames(mat)[order(-m)][seq_len(min(top, nrow(mat)))]
}

# Feature levels, bottom to top: by mean value (largest on top) or by the
# order of a hierarchical clustering of the feature profiles.
#' @importFrom stats aggregate dist hclust na.pass
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

# Tile centre and width of each time within each panel.
.add_tile_edges <- function(df, visit) {
    edges <- lapply(split(df$.time, df$.panel), function(t) {
        t <- sort(unique(t))
        e <- if (visit) list(lower = t - 0.5, upper = t + 0.5) else
            .tile_edges(t)
        data.frame(.time = t, .lower = e$lower, .upper = e$upper)
    })
    edges <- do.call(rbind, Map(cbind, .panel = names(edges), edges))
    df <- merge(df, edges, by = c(".panel", ".time"), sort = FALSE)
    df$.mid <- (df$.lower + df$.upper) / 2
    df$.width <- df$.upper - df$.lower
    df
}

#' @importFrom ggplot2 ggplot aes .data geom_tile labs facet_wrap vars
#'   scale_fill_gradient scale_fill_gradient2 scale_x_continuous
.plot_tiles <- function(df, time.col, assay.type, value, level, panel.col,
        levels, ...) {
    p <- ggplot(df, aes(x = .data$.mid, y = .data$.feature)) +
        geom_tile(aes(width = .data$.width, fill = .data$.value), ...)
    lab <- switch(value, raw = assay.type,
        change = paste(assay.type, "change"),
        occupancy = if (level == "group") "occupancy" else "present")
    p <- p + if (value == "change") {
        scale_fill_gradient2(low = "#b2182b", mid = "grey95",
            high = "#2166ac", na.value = "white")
    } else {
        scale_fill_gradient(low = "grey95", high = "#2c6fbb",
            na.value = "white", limits = if (value == "occupancy") c(0, 1))
    }
    if (!is.null(levels)) {
        p <- p + scale_x_continuous(breaks = seq_along(levels),
            labels = levels)
    }
    if (!is.null(panel.col)) {
        p <- p + facet_wrap(vars(.data$.panel))
    }
    p + labs(x = time.col, y = NULL, fill = lab)
}
