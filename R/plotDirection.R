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
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param dimred \code{Character scalar}. Specifies the
#'   \code{reducedDim(x)} to use.
#'
#' @param ncomponents \code{Numeric} or \code{Character vector}. The two
#'   components of \code{dimred}, by index or name. (Default: \code{c(1, 2)})
#'
#' @param bins \code{Integer scalar}. Number of direction sectors.
#'   (Default: \code{16})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the moves into panels. Each move takes
#'   the value of its earlier sample. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_histogram]{geom_histogram}}.
#'
#' @details
#' Angles are measured counterclockwise from the first component (east),
#' so 90 degrees points along the second component. Moves of zero length
#' have no direction and are left out.
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
#'     dimred = "PCoA", group.by = "StudyIdentifier")
NULL

#' @rdname plotDirection
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom methods setMethod
#' @export
setMethod("plotDirection", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, dimred, ncomponents = c(1, 2), bins = 16,
            group.by = NULL, ...) {
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        .time_type(x, time.col)
        .check_column(x, group.by)
        if (!(is.numeric(bins) && length(bins) == 1L && bins >= 2)) {
            stop("'bins' must be a single number of at least 2.",
                call. = FALSE)
        }
        labels <- paste0(dimred, ":", ncomponents)
        hits <- lapply(labels, .match_dimred_component, x = x)
        if (length(labels) != 2L ||
                any(vapply(hits, is.null, logical(1)))) {
            stop("'dimred' must be in reducedDimNames(x) and ",
                "'ncomponents' must give two of its components.",
                call. = FALSE)
        }
        rd <- reducedDim(x, dimred)
        cols <- unique(c(time.col, group, group.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df$.x <- as.numeric(rd[, hits[[1L]]$j])
        df$.y <- as.numeric(rd[, hits[[2L]]$j])
        df <- .order_time(df, time.col, group)
        n <- nrow(df)
        i <- which(df[[group]][-1L] == df[[group]][-n])
        dx <- df$.x[i + 1L] - df$.x[i]
        dy <- df$.y[i + 1L] - df$.y[i]
        moves <- data.frame(
            .angle = (atan2(dy, dx) * 180 / pi) %% 360,
            .length = sqrt(dx^2 + dy^2))
        if (!is.null(group.by)) {
            moves[[group.by]] <- df[[group.by]][i]
        }
        moves <- moves[moves$.length > 0, , drop = FALSE]
        if (nrow(moves) == 0L) {
            stop("No subject moves between consecutive samples.",
                call. = FALSE)
        }
        .plot_direction(moves, labels, as.integer(bins), group.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_histogram coord_polar labs
#'   scale_x_continuous facet_wrap vars
.plot_direction <- function(df, labels, bins, group.by, ...) {
    # Sectors are centred on their directions: the first spans east.
    w <- 360 / bins
    df$.angle[df$.angle >= 360 - w / 2] <- df$.angle[df$.angle >=
        360 - w / 2] - 360
    p <- ggplot(df, aes(x = .data$.angle)) +
        geom_histogram(breaks = seq(-w / 2, 360 - w / 2, by = w),
            fill = "#2c6fbb", colour = "white", ...) +
        scale_x_continuous(limits = c(-w / 2, 360 - w / 2),
            breaks = c(0, 90, 180, 270),
            labels = c(labels[1L], labels[2L], paste0("-", labels[1L]),
                paste0("-", labels[2L]))) +
        coord_polar(start = -pi / 2 - w / 2 * pi / 180, direction = -1) +
        labs(x = NULL, y = "moves")
    if (!is.null(group.by)) {
        p <- p + facet_wrap(vars(.data[[group.by]]))
    }
    p
}
