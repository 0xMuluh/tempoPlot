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
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
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
#' @param top \code{Integer scalar}. Number of features with the largest
#'   mean value in \code{assay.type} to draw. (Default: \code{5})
#'
#' @param coord \code{Character scalar}. Draws time on a straight
#'   (\code{"cartesian"}) or circular (\code{"polar"}) axis.
#'   (Default: \code{"cartesian"})
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
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom methods setMethod
#' @export
setMethod("plotRanks", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts", top = 5,
            coord = c("cartesian", "polar"), ...) {
        coord <- match.arg(coord)
        .check_longitudinal(x, time.col, group)
        .time_type(x, time.col)
        if (!.is_one_name(assay.type) || !assay.type %in% assayNames(x)) {
            stop("'assay.type' must be one of assayNames(x).",
                call. = FALSE)
        }
        if (!(is.numeric(top) && length(top) == 1L && top >= 1)) {
            stop("'top' must be a single number of at least 1.",
                call. = FALSE)
        }
        mat <- as.matrix(assay(x, assay.type))
        # matrix() keeps the dimensions when one feature is left.
        ranks <- matrix(apply(-mat, 2L, base::rank, ties.method = "average"),
            nrow(mat), dimnames = dimnames(mat))
        keep <- order(rowMeans(mat, na.rm = TRUE),
            decreasing = TRUE)[seq_len(min(top, nrow(mat)))]
        ranks <- ranks[keep, , drop = FALSE]

        time <- colData(x)[[time.col]]
        subject <- if (is.null(group)) rep("", ncol(x)) else
            as.character(colData(x)[[group]])
        df <- data.frame(
            .feature = factor(rep(rownames(ranks), times = ncol(x)),
                levels = rownames(ranks)),
            .panel = rep(subject, each = nrow(ranks)),
            .time = rep(time, each = nrow(ranks)),
            .rank = as.vector(ranks))
        df <- df[order(df$.feature, df$.panel, .time_numeric(df$.time)), ,
            drop = FALSE]
        .plot_ranks(df, time.col, coord, !is.null(group), ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs
#'   scale_y_reverse coord_polar facet_wrap vars
.plot_ranks <- function(df, time.col, coord, panels, ...) {
    p <- ggplot(df, aes(x = .data$.time, y = .data$.rank,
            colour = .data$.feature,
            group = interaction(.data$.feature, .data$.panel))) +
        geom_line(...) +
        geom_point(size = 1) +
        scale_y_reverse() +
        labs(x = time.col, y = "rank", colour = "feature")
    if (coord == "polar") {
        p <- p + coord_polar()
    }
    if (panels) {
        p <- p + facet_wrap(vars(.data$.panel))
    }
    p
}
