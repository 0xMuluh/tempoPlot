#' Plot features over time as horizon strips
#'
#' \code{plotHorizon()} draws one strip per feature. Each value is centred,
#' cut into \code{bands} bands of equal height and the bands are overlaid,
#' so that deviations above and below the centre show as darker blues and
#' reds in a short strip.
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
#' @param bands \code{Integer scalar}. Number of bands on each side of the
#'   centre. (Default: \code{3})
#'
#' @param centre \code{Character scalar}. Centres each series on its median
#'   (\code{"median"}), its first value (\code{"baseline"}) or its mean
#'   (\code{"mean"}). (Default: \code{"median"})
#'
#' @param level \code{Character scalar}. Draws one column of strips per
#'   subject (\code{"subject"}) or averages the subjects of each
#'   \code{group.by} group at each time (\code{"group"}).
#'   (Default: \code{"subject"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that groups subjects when \code{level = "group"}. If
#'   \code{NULL}, all subjects form one group. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_ribbon]{geom_ribbon}}.
#'
#' @details
#' The band height of a feature is its largest absolute deviation from the
#' centre, over all subjects, divided by \code{bands}; strips of a feature
#' are comparable across subjects, not across features. Visits are evenly
#' spaced.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Fink I, et al. (2022). BiomeHorizon: visualizing microbiome time series
#' data in R. \emph{mSystems} 7(3):e01380-21.
#'
#' @name plotHorizon
#'
#' @examples
#' # 16 species in three bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "relabundance")
#' plotHorizon(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     assay.type = "relabundance", top = 8)
#'
#' # Deviation from the first sample, averaged over the bioreactors
#' plotHorizon(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     assay.type = "relabundance", top = 8, centre = "baseline",
#'     level = "group")
NULL

#' @rdname plotHorizon
#' @importFrom SummarizedExperiment colData
#' @importFrom stats ave median
#' @importFrom methods setMethod
#' @export
setMethod("plotHorizon", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            features = NULL, top = 20, bands = 3,
            centre = c("median", "baseline", "mean"),
            level = c("subject", "group"), group.by = NULL, ...) {
        centre <- match.arg(centre)
        level <- match.arg(level)
        if (!(is.numeric(bands) && length(bands) == 1L && bands >= 1)) {
            stop("'bands' must be a single number of at least 1.",
                call. = FALSE)
        }
        bands <- as.integer(bands)
        df <- .feature_table(x, time.col, group, assay.type, features, top,
            "raw", level, group.by)
        type <- .time_type(x, time.col)
        # Features in selection order, the most abundant at the top.
        df$.feature <- factor(df$.feature, levels = unique(df$.feature))
        df <- df[order(df$.feature, df$.panel, df$.time), , drop = FALSE]

        series <- interaction(df$.feature, df$.panel, drop = TRUE)
        mid <- switch(centre,
            median = ave(df$.value, series, FUN = function(v)
                median(v, na.rm = TRUE)),
            mean = ave(df$.value, series, FUN = function(v)
                mean(v, na.rm = TRUE)),
            baseline = ave(df$.value, series, FUN = function(v) v[1L]))
        df$.value <- df$.value - mid
        df <- .horizon_bands(df, bands)
        if (type == "date") {
            df$.time <- as.Date(df$.time, origin = "1970-01-01")
        }
        levels <- if (type == "visit") levels(colData(x)[[time.col]])
        panel <- !is.null(if (level == "subject") group else group.by)
        .plot_horizon(df, time.col, bands, panel, levels, ...)
    }
)

# One row per value and band: .band is +k above and -k below the centre,
# .y the part of the deviation inside band k, in band heights (0 to 1).
.horizon_bands <- function(df, bands) {
    height <- ave(abs(df$.value), df$.feature, FUN = function(v) {
        m <- suppressWarnings(max(v, na.rm = TRUE))
        if (is.finite(m) && m > 0) m / bands else 1
    })
    z <- df$.value / height
    res <- lapply(seq_len(bands), function(k) {
        up <- pmin(pmax(z - (k - 1L), 0), 1)
        down <- pmin(pmax(-z - (k - 1L), 0), 1)
        rbind(cbind(df, .band = k, .y = up),
            cbind(df, .band = -k, .y = down))
    })
    res <- do.call(rbind, res)
    # Legend from the highest band above the centre to the lowest below.
    res$.band <- factor(res$.band, levels = c(rev(seq_len(bands)),
        -seq_len(bands)))
    res
}

#' @importFrom ggplot2 ggplot aes .data geom_ribbon labs facet_grid vars
#'   scale_fill_manual scale_x_continuous theme element_blank element_text
#' @importFrom grDevices colorRampPalette
.plot_horizon <- function(df, time.col, bands, panel, levels, ...) {
    colours <- c(
        colorRampPalette(c("#f4c7c3", "#b2182b"))(bands),
        colorRampPalette(c("#c6dbef", "#2166ac"))(bands))
    names(colours) <- c(-seq_len(bands), seq_len(bands))
    p <- ggplot(df, aes(x = .data$.time, ymin = 0, ymax = .data$.y,
            fill = .data$.band,
            group = interaction(.data$.band, .data$.panel))) +
        geom_ribbon(...) +
        scale_fill_manual(values = colours, drop = FALSE) +
        facet_grid(rows = vars(.data$.feature),
            cols = if (panel) vars(.data$.panel), switch = "y") +
        labs(x = time.col, y = NULL, fill = "band") +
        theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 1))
    if (!is.null(levels)) {
        p <- p + scale_x_continuous(breaks = seq_along(levels),
            labels = levels)
    }
    p
}
