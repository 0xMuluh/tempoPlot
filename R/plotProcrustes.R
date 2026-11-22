#' Plot two ordinations superimposed by Procrustes rotation
#'
#' \code{plotProcrustes()} rotates a second ordination of the same samples
#' onto the first and draws both, with a segment from each sample in the
#' first to the same sample in the second. Short segments mean the two
#' ordinations agree.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object, such as a \code{TreeSummarizedExperiment}, or a
#'   \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}.
#'
#' @param dimred \code{Character scalar}. Specifies the first
#'   \code{reducedDim}.
#'
#' @param dimred2 \code{Character scalar}. Specifies the second
#'   \code{reducedDim}, of \code{x} or of \code{experiment2}. For a
#'   \code{MultiAssayExperiment}, it defaults to \code{dimred}.
#'
#' @param ncomponents \code{Numeric scalar} or \code{vector}. The
#'   components of each ordination used in the rotation: the first
#'   \code{ncomponents}, or the components given by index or name. The
#'   first two are drawn. (Default: \code{2})
#'
#' @param colour.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that colours the points. (Default: \code{NULL})
#'
#' @param experiment,experiment2 \code{Character scalar}. For a
#'   \code{MultiAssayExperiment}, the experiments holding \code{dimred} and
#'   \code{dimred2}. Samples are matched by their primary names.
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Both ordinations are centred and scaled to unit sum of squares, and the
#' second is rotated to the first (symmetric Procrustes analysis). The
#' subtitle gives the Procrustes statistic \eqn{m^2}, the sum of squared
#' distances between matched samples, from 0 (same configuration) to 1.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Peres-Neto PR, Jackson DA (2001). How well do multivariate data sets
#' match? The advantages of a Procrustean superimposition approach over
#' the Mantel test. \emph{Oecologia} 129(2):169-178.
#'
#' @name plotProcrustes
#'
#' @examples
#' # PCoA of Bray-Curtis against PCA of clr abundances
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "relabundance")
#' tse <- mia::transformAssay(tse, method = "clr", pseudocount = 1,
#'     assay.type = "counts")
#' tse <- mia::addMDS(tse, assay.type = "relabundance", method = "bray",
#'     name = "PCoA")
#' tse <- scater::runPCA(tse, assay.type = "clr", ncomponents = 5,
#'     name = "PCA")
#' plotProcrustes(tse, dimred = "PCoA", dimred2 = "PCA",
#'     colour.by = "StudyIdentifier")
NULL

#' @rdname plotProcrustes
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotProcrustes", signature = c(x = "SummarizedExperiment"),
    function(x, dimred, dimred2, ncomponents = 2, colour.by = NULL, ...) {
        .check_column(x, colour.by)
        a <- .get_reduced_dim(x, dimred, ncomponents, "dimred")
        b <- .get_reduced_dim(x, dimred2, ncomponents, "dimred2")
        colour <- if (!is.null(colour.by)) colData(x)[[colour.by]]
        .plot_procrustes(a, b, dimred, dimred2, colour, colour.by, ...)
    }
)

#' @rdname plotProcrustes
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom methods setMethod
#' @export
setMethod("plotProcrustes", signature = c(x = "MultiAssayExperiment"),
    function(x, experiment, experiment2, dimred, dimred2 = dimred,
            ncomponents = 2, colour.by = NULL, ...) {
        exps <- names(experiments(x))
        for (arg in c("experiment", "experiment2")) {
            if (!.is_one_name(get(arg)) || !get(arg) %in% exps) {
                stop("'", arg, "' must be one experiment of 'x'.",
                    call. = FALSE)
            }
        }
        cd <- colData(x)
        if (!is.null(colour.by) && !colour.by %in% names(cd)) {
            stop("'colour.by' must be one column name of colData(x).",
                call. = FALSE)
        }
        e1 <- x[[experiment]]
        e2 <- x[[experiment2]]
        a <- .get_reduced_dim(e1, dimred, ncomponents, "dimred")
        b <- .get_reduced_dim(e2, dimred2, ncomponents, "dimred2")
        map <- sampleMap(x)
        primary <- function(e, cols) {
            m <- map[map$assay == e, , drop = FALSE]
            m$primary[match(cols, m$colname)]
        }
        p1 <- primary(experiment, colnames(e1))
        p2 <- primary(experiment2, colnames(e2))
        common <- intersect(p1, p2)
        if (length(common) < 3L) {
            stop("The experiments share fewer than three samples.",
                call. = FALSE)
        }
        a <- a[match(common, p1), , drop = FALSE]
        b <- b[match(common, p2), , drop = FALSE]
        colour <- if (!is.null(colour.by)) cd[common, colour.by]
        .plot_procrustes(a, b, paste(experiment, dimred),
            paste(experiment2, dimred2), colour, colour.by, ...)
    }
)

# Symmetric Procrustes: centre, scale to unit sum of squares and rotate b
# onto a. Returns both configurations and m2.
.procrustes <- function(a, b) {
    std <- function(m) {
        m <- scale(m, scale = FALSE)
        m / sqrt(sum(m^2))
    }
    a <- std(a)
    b <- std(b)
    s <- svd(crossprod(a, b))
    rot <- s$v %*% t(s$u)
    list(a = a, b = b %*% rot, m2 = 1 - sum(s$d)^2)
}

#' @importFrom ggplot2 ggplot aes .data geom_point geom_segment labs arrow
#'   unit coord_equal
.plot_procrustes <- function(a, b, lab1, lab2, colour, colour.by, ...) {
    if (ncol(a) < 2L || ncol(b) < 2L) {
        stop("'ncomponents' must select at least two components.",
            call. = FALSE)
    }
    fit <- .procrustes(a, b)
    df <- data.frame(.x = fit$a[, 1L], .y = fit$a[, 2L],
        .x2 = fit$b[, 1L], .y2 = fit$b[, 2L])
    if (!is.null(colour.by)) {
        df[[colour.by]] <- colour
    }
    ggplot(df, aes(colour = if (!is.null(colour.by)) .data[[colour.by]])) +
        geom_segment(aes(x = .data$.x, y = .data$.y, xend = .data$.x2,
            yend = .data$.y2), arrow = arrow(length = unit(0.12, "cm")),
            alpha = 0.6) +
        geom_point(aes(x = .data$.x, y = .data$.y), ...) +
        .scale_by(df, colour.by, "colour") +
        coord_equal() +
        labs(x = "dimension 1", y = "dimension 2", colour = colour.by,
            subtitle = sprintf("%s (points) to %s (arrows); m2 = %.3f",
                lab1, lab2, fit$m2)) +
        .theme_tempoPlot()
}
