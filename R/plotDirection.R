#' Plot directions of movement in an ordination
#'
#' \code{plotDirection()} draws a rose diagram of the directions in which
#' subjects move between consecutive samples in two components of a
#' reduced dimension.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object, such as a \code{TreeSummarizedExperiment}.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits. Time differences, as in
#'   \code{max.gap}, are in days for \code{Date} and hours for
#'   \code{POSIXct} times.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param dimred \code{Character scalar}. Specifies the
#'   \code{reducedDim(x)} to use.
#'
#' @param ncomponents \code{Numeric scalar} or \code{vector}. The two
#'   components of \code{dimred}: \code{2} for the first two, or two
#'   components by index or name. (Default: \code{2})
#'
#' @param bins \code{Integer scalar}. Number of direction sectors.
#'   (Default: \code{16})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the moves into panels. Each move takes
#'   the value of its later sample. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Consecutive samples of a subject
#'   further apart than this, in \code{time.col} units, are not a move.
#'   For visits, \code{NULL} skips moves across visits without a sample;
#'   otherwise \code{NULL} counts all. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_histogram]{geom_histogram}}.
#'
#' @details
#' Angles are measured counterclockwise from the first component (east),
#' so 90 degrees points along the second component. Moves of zero length
#' have no direction and are left out. The radius is on a square-root
#' scale, so that the area of a sector follows its number of moves.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotDirection
#'
#' @examples
#' # Directions in a PCoA of three bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "relabundance")
#' tse <- mia::addMDS(tse, assay.type = "relabundance", method = "bray",
#'     name = "PCoA")
#' plotDirection(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     dimred = "PCoA", facet.by = "StudyIdentifier")
NULL

#' @rdname plotDirection
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotDirection", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, dimred, ncomponents = 2, bins = 16,
            facet.by = NULL, max.gap = NULL, ...) {
        .require_column(x, group)
        .check_longitudinal(x, time.col, group)
        .check_column(x, facet.by)
        .check_number(bins, 2)
        rd <- .get_reduced_dim(x, dimred, ncomponents, n = 2L)
        cols <- unique(c(time.col, group, facet.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df$.x <- rd[, 1L]
        df$.y <- rd[, 2L]
        df <- .order_time(df, time.col, group)
        df <- .break_gaps(df, time.col, group, max.gap)
        i <- .next_in_subject(df[[group]], segment = df$.segment)
        dx <- df$.x[i + 1L] - df$.x[i]
        dy <- df$.y[i + 1L] - df$.y[i]
        moves <- data.frame(
            .angle = (atan2(dy, dx) * 180 / pi) %% 360,
            .length = sqrt(dx^2 + dy^2))
        if (!is.null(facet.by)) {
            moves[[facet.by]] <- df[[facet.by]][i + 1L]
        }
        moves <- moves[moves$.length > 0, , drop = FALSE]
        if (nrow(moves) == 0L) {
            stop("No subject moves between consecutive samples.",
                call. = FALSE)
        }
        .plot_direction(moves, colnames(rd), as.integer(bins), facet.by,
            ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_histogram coord_polar labs
#'   scale_x_continuous scale_y_sqrt facet_wrap vars
#' @importFrom utils modifyList
.plot_direction <- function(df, labels, bins, facet.by, ...) {
    # Sectors are centred on their directions: the first spans east.
    w <- 360 / bins
    df$.angle[df$.angle >= 360 - w / 2] <- df$.angle[df$.angle >=
        360 - w / 2] - 360
    p <- ggplot(df, aes(x = .data$.angle)) +
        do.call(geom_histogram, modifyList(list(
            breaks = seq(-w / 2, 360 - w / 2, by = w), fill = .okabe_ito[[1L]],
            colour = "white"), list(...))) +
        scale_x_continuous(limits = c(-w / 2, 360 - w / 2),
            breaks = c(0, 90, 180, 270),
            labels = c(labels[1L], labels[2L], paste0("-", labels[1L]),
                paste0("-", labels[2L]))) +
        scale_y_sqrt() +
        coord_polar(start = -pi / 2 - w / 2 * pi / 180, direction = -1,
            clip = "off") +
        labs(x = NULL, y = "moves")
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    p + .theme_tempoPlot()
}
