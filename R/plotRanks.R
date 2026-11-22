#' Plot the ranks of features over time
#'
#' \code{plotRanks()} draws the abundance rank of the most abundant
#' features in each sample as bump lines, one panel per subject. In polar
#' coordinates the plot becomes a rank clock.
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
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of
#'   abundances. To rank taxa such as phyla, agglomerate \code{x} first with
#'   \code{\link[mia:agglomerate-methods]{mia::agglomerateByRank()}}.
#'   (Default: \code{"counts"})
#'
#' @param features \code{Character vector}. Features to draw. If
#'   \code{NULL}, draws the \code{top} features with the best mean rank.
#'   (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features to draw when
#'   \code{features} is \code{NULL}. (Default: \code{5})
#'
#' @param coord \code{Character scalar}. Draws time on a straight
#'   (\code{"cartesian"}) or circular (\code{"polar"}) axis.
#'   (Default: \code{"cartesian"})
#'
#' @param max.gap \code{Numeric scalar}. Breaks a line where consecutive
#'   samples of a subject are further apart than this, in \code{time.col}
#'   units. For visits, \code{NULL} breaks lines at visits without a
#'   sample; otherwise \code{NULL} never breaks lines. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' Ranks are among all features of a sample, 1 the most abundant; ties get
#' their mean rank.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Collins SL, et al. (2008). Rank clocks and plant community dynamics.
#' \emph{Ecology} 89(12):3534-3541.
#'
#' @name plotRanks
#'
#' @examples
#' # 16 species in three bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- minimalgut
#' plotRanks(tse, time.col = "Time.hr", group = "StudyIdentifier")
#'
#' # As rank clocks
#' plotRanks(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     coord = "polar")
#'
#' # Phyla of the subjects sampled five times
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006[, !duplicated(t(SummarizedExperiment::assay(
#'     hitchip1006, "counts")))]
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]
#' tse <- mia::agglomerateByRank(tse, rank = "Phylum")
#' plotRanks(tse, time.col = "time", group = "subject", top = 3)
NULL

#' @rdname plotRanks
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotRanks", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            features = NULL, top = 5, coord = c("cartesian", "polar"),
            max.gap = NULL, ...) {
        coord <- match.arg(coord)
        .check_longitudinal(x, time.col, group)
        mat <- .get_assay(x, assay.type)
        # matrix() keeps the dimensions when one feature is left.
        ranks <- matrix(apply(-mat, 2L, base::rank, ties.method = "average"),
            nrow(mat), dimnames = dimnames(mat))
        ranks <- ranks[.select_features(-ranks, features, top), ,
            drop = FALSE]

        cd <- data.frame(.panel = if (is.null(group)) "" else
            as.character(colData(x)[[group]]),
            .time = colData(x)[[time.col]], .i = seq_len(ncol(x)))
        cd <- .order_time(cd, ".time", ".panel")
        cd <- .break_gaps(cd, ".time", ".panel", max.gap)
        df <- data.frame(
            .feature = factor(rep(rownames(ranks), times = ncol(x)),
                levels = rownames(ranks)),
            cd[rep(seq_len(ncol(x)), each = nrow(ranks)), c(".panel",
                ".time", ".segment")],
            .rank = as.vector(ranks[, cd$.i, drop = FALSE]))
        .plot_ranks(df, time.col, coord, !is.null(group), ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs
#'   scale_y_reverse scale_x_discrete coord_polar facet_wrap vars
.plot_ranks <- function(df, time.col, coord, panels, ...) {
    p <- ggplot(df, aes(x = .data$.time, y = .data$.rank,
            colour = .data$.feature,
            group = interaction(.data$.feature, .data$.panel,
                .data$.segment))) +
        geom_line(...) +
        geom_point(size = 1) +
        scale_y_reverse() +
        .scale_categorical("colour") +
        labs(x = time.col, y = "rank", colour = "feature") +
        .theme_tempoPlot()
    if (is.factor(df$.time)) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    if (coord == "polar") {
        p <- p + coord_polar()
    }
    if (panels) {
        p <- p + facet_wrap(vars(.data$.panel))
    }
    p
}
