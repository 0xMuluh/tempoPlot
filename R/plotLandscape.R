#' Plot the density landscape of a variable
#'
#' \code{plotLandscape()} draws the density of a variable over all
#' samples, or its potential, the negative log density. Peaks of the
#' density, valleys of the potential, suggest alternative states.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param y \code{Character scalar}. Specifies the variable: a column of
#'   \code{colData(x)}, a row of \code{assay(x, assay.type)}, or
#'   \code{"<dimred>:<component>"}.
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to
#'   retrieve \code{y} from when \code{y} is a feature name.
#'   (Default: \code{NULL})
#'
#' @param scale \code{Character scalar}. Draws the density
#'   (\code{"density"}) or the potential (\code{"potential"}).
#'   (Default: \code{"density"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the samples into curves.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[stats:density]{density}}, such as \code{bw} or
#'   \code{adjust}.
#'
#' @details
#' Densities are Gaussian kernel estimates; the bandwidth sets how many
#' states can appear. Each sample counts once, so subjects with more
#' samples weigh more. Ticks on the x axis mark the samples.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Lahti L, et al. (2014). Tipping elements in the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.
#'
#' @name plotLandscape
#'
#' @examples
#' # Prevotella abundance at the first sample of each subject
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006[, hitchip1006$time == 0]
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' tse <- mia::transformAssay(tse, assay.type = "relabundance",
#'     method = "log10", pseudocount = 1e-4, name = "log10")
#' plotLandscape(tse, y = "Prevotella melaninogenica et rel.",
#'     assay.type = "log10", scale = "potential")
NULL

#' @rdname plotLandscape
#' @importFrom SummarizedExperiment colData
#' @importFrom stats density
#' @importFrom methods setMethod
#' @export
setMethod("plotLandscape", signature = c(x = "SummarizedExperiment"),
    function(x, y, assay.type = NULL, scale = c("density", "potential"),
            group.by = NULL, ...) {
        scale <- match.arg(scale)
        .check_column(x, group.by)
        value <- unname(.resolve_y(x, y, assay.type))
        by <- if (is.null(group.by)) rep("", ncol(x)) else
            as.character(colData(x)[[group.by]])
        keep <- !is.na(value) & !is.na(by)
        groups <- split(value[keep], by[keep])
        groups <- groups[lengths(groups) >= 2L]
        if (length(groups) == 0L) {
            stop("'y' needs at least two values in a group.",
                call. = FALSE)
        }
        curves <- do.call(rbind, lapply(seq_along(groups), function(i) {
            d <- density(groups[[i]], ...)
            data.frame(.by = names(groups)[i], .x = d$x, .y = d$y)
        }))
        if (scale == "potential") {
            curves <- curves[curves$.y > 0, , drop = FALSE]
            curves$.y <- -log(curves$.y)
        }
        samples <- data.frame(.by = by[keep], .x = value[keep])
        .plot_landscape(curves, samples, y, scale, group.by)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_rug labs
.plot_landscape <- function(curves, samples, y, scale, group.by) {
    colour <- if (!is.null(group.by)) aes(colour = .data$.by)
    ggplot(curves, aes(x = .data$.x, y = .data$.y, group = .data$.by)) +
        geom_line(colour, linewidth = 1) +
        geom_rug(aes(x = .data$.x, colour = if (!is.null(group.by))
            .data$.by), data = samples, inherit.aes = FALSE, sides = "b",
            alpha = 0.3) +
        labs(x = y, y = if (scale == "density") "density" else
            "potential (-log density)", colour = group.by)
}
