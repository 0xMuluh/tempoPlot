#' Plot trajectories of subjects in two dimensions
#'
#' \code{plotPath()} joins each subject's samples in time order in a
#' two-dimensional space: two components of an ordination, or two variables
#' such as features or colData columns.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param dimred \code{Character scalar}. Specifies the
#'   \code{reducedDim(x)} to plot. (Default: \code{NULL})
#'
#' @param ncomponents \code{Numeric} or \code{Character vector}. The two
#'   components of \code{dimred} to plot, by index or name.
#'   (Default: \code{c(1, 2)})
#'
#' @param x.var,y.var \code{Character scalar}. Variables to plot instead of
#'   \code{dimred}: a column of \code{colData(x)} or a feature of
#'   \code{assay.type}. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay that holds
#'   features named in \code{x.var} or \code{y.var}. (Default: \code{NULL})
#'
#' @param origin \code{Character scalar}. Whether paths start where the
#'   first sample lies (\code{"none"}) or at zero, by subtracting each
#'   subject's first sample (\code{"baseline"}). (Default: \code{"none"})
#'
#' @param arrows \code{Logical scalar}. Whether paths end in an arrow
#'   pointing forward in time. (Default: \code{TRUE})
#'
#' @param colour.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that colours the paths. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Breaks a path where consecutive
#'   samples of a subject are further apart than this, in \code{time.col}
#'   units. \code{NULL} never breaks paths. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_path]{geom_path}}.
#'
#' @details
#' Give either \code{dimred} or both \code{x.var} and \code{y.var}.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotPath
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
#' # Subjects sampled five times
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]
#'
#' # Paths of subjects in a PCoA of Bray-Curtis dissimilarities
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' tse <- mia::addMDS(tse, assay.type = "relabundance", method = "bray",
#'     name = "PCoA")
#' plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
#'     colour.by = "subject")
#'
#' # The same paths, all starting at zero, to compare their directions
#' plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
#'     origin = "baseline", colour.by = "subject")
#'
#' # Two genera against each other
#' plotPath(tse, time.col = "time", group = "subject",
#'     x.var = "Bacteroides vulgatus et rel.",
#'     y.var = "Prevotella melaninogenica et rel.",
#'     assay.type = "relabundance", colour.by = "subject")
NULL

#' @rdname plotPath
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom methods setMethod
#' @export
setMethod("plotPath", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, dimred = NULL, ncomponents = c(1, 2),
            x.var = NULL, y.var = NULL, assay.type = NULL,
            origin = c("none", "baseline"), arrows = TRUE, colour.by = NULL,
            max.gap = NULL, ...) {
        origin <- match.arg(origin)
        if (!(isTRUE(arrows) || isFALSE(arrows))) {
            stop("'arrows' must be TRUE or FALSE.", call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        # Errors on character and other unsupported time types.
        .time_type(x, time.col)
        .check_column(x, colour.by)
        if (is.null(dimred) == (is.null(x.var) || is.null(y.var))) {
            stop("Give either 'dimred' or both 'x.var' and 'y.var'.",
                call. = FALSE)
        }

        cols <- unique(c(time.col, group, colour.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        if (is.null(dimred)) {
            df$.x <- unname(.resolve_y(x, x.var, assay.type, "x.var"))
            df$.y <- unname(.resolve_y(x, y.var, assay.type, "y.var"))
            labels <- c(x.var, y.var)
        } else {
            labels <- paste0(dimred, ":", ncomponents)
            hits <- lapply(labels, .match_dimred_component, x = x)
            if (length(labels) != 2L || any(vapply(hits, is.null,
                    logical(1)))) {
                stop("'dimred' must be in reducedDimNames(x) and ",
                    "'ncomponents' must give two of its components.",
                    call. = FALSE)
            }
            rd <- reducedDim(x, dimred)
            df$.x <- as.numeric(rd[, hits[[1L]]$j])
            df$.y <- as.numeric(rd[, hits[[2L]]$j])
        }
        df <- .order_time(df, time.col, group)
        df <- .break_gaps(df, time.col, group, max.gap)
        if (origin == "baseline") {
            first <- function(v) ave(v, df[[group]], FUN = function(u) u[1L])
            df$.x <- df$.x - first(df$.x)
            df$.y <- df$.y - first(df$.y)
        }
        df$.line <- interaction(df[[group]], df$.segment, drop = TRUE)

        .plot_path(df, labels, arrows, colour.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_path geom_point labs arrow unit
.plot_path <- function(df, labels, arrows, colour.by, ...) {
    ggplot(df, aes(x = .data$.x, y = .data$.y,
            colour = if (!is.null(colour.by)) .data[[colour.by]])) +
        geom_path(aes(group = .data$.line),
            arrow = if (arrows) arrow(length = unit(0.2, "cm")), ...) +
        geom_point(size = 1) +
        labs(x = labels[1L], y = labels[2L], colour = colour.by)
}
