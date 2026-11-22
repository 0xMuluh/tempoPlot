#' Plot dissimilarity against overlap of sample pairs
#'
#' \code{plotDOC()} draws the dissimilarity of each pair of samples
#' against their overlap, with a smooth: the dissimilarity-overlap curve
#' (DOC). A falling curve at high overlap suggests that the same taxa
#' behave alike in different communities.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param overlap,dissimilarity \code{Character scalar}. Specify the
#'   sample by sample matrices in \code{metadata(x)}, as stored by
#'   \code{\link[mia:getDissimilarity]{mia::addDissimilarity()}}: the
#'   overlap on the x axis and the dissimilarity on the y axis.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects, needed for
#'   \code{pairs = "within"} and \code{"between"}. (Default: \code{NULL})
#'
#' @param pairs \code{Character scalar}. Draws all pairs of samples
#'   (\code{"all"}), pairs from the same subject (\code{"within"}) or from
#'   different subjects (\code{"between"}). (Default: \code{"all"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Each pair of samples is drawn once, with the robust LOWESS curve of
#' Bashan et al. (2016); above 1000 pairs, points are drawn with
#' transparency. In Bashan et al. (2016), overlap is mia's
#' \code{"overlap"} of relative abundances and dissimilarity is the root
#' Jensen-Shannon divergence of the shared taxa, renormalised; mia's
#' \code{"jsd"} uses all taxa.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Bashan A, et al. (2016). Universality of human microbial dynamics.
#' \emph{Nature} 534(7606):259-262.
#'
#' @name plotDOC
#'
#' @examples
#' # Samples of four artificial guts
#' data("SilvermanAGutData", package = "miaTime")
#' tse <- SilvermanAGutData
#' # Five samples have no reads
#' tse <- tse[, colSums(SummarizedExperiment::assay(tse)) > 0]
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' tse <- mia::addDissimilarity(tse, method = "overlap",
#'     assay.type = "relabundance")
#' tse <- mia::addDissimilarity(tse, method = "jsd",
#'     assay.type = "relabundance")
#' plotDOC(tse, overlap = "overlap", dissimilarity = "jsd",
#'     group = "Vessel", pairs = "between")
NULL

#' @rdname plotDOC
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotDOC", signature = c(x = "SummarizedExperiment"),
    function(x, overlap, dissimilarity, group = NULL,
            pairs = c("all", "within", "between"), ...) {
        pairs <- match.arg(pairs)
        .check_column(x, group)
        if (pairs != "all" && is.null(group)) {
            stop("pairs = \"", pairs, "\" needs 'group'.", call. = FALSE)
        }
        o <- .get_dissimilarity(x, overlap, "overlap")
        d <- .get_dissimilarity(x, dissimilarity)
        ij <- which(lower.tri(o), arr.ind = TRUE)
        if (pairs != "all") {
            subject <- as.character(colData(x)[[group]])
            same <- subject[ij[, 1L]] == subject[ij[, 2L]]
            ij <- ij[if (pairs == "within") same else !same, ,
                drop = FALSE]
        }
        df <- data.frame(.overlap = o[ij], .dissimilarity = d[ij])
        df <- df[!is.na(df$.overlap) & !is.na(df$.dissimilarity), ,
            drop = FALSE]
        if (nrow(df) == 0L) {
            stop("No pairs of samples to draw.", call. = FALSE)
        }
        .plot_doc(df, overlap, dissimilarity, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_point geom_line labs
#' @importFrom stats lowess
#' @importFrom utils modifyList
.plot_doc <- function(df, overlap, dissimilarity, ...) {
    n <- nrow(df)
    alpha <- max(0.05, min(1, 1000 / n))
    p <- ggplot(df, aes(x = .data$.overlap, y = .data$.dissimilarity)) +
        do.call(geom_point, modifyList(list(alpha = alpha), list(...)))
    if (n >= 10L) {
        curve <- as.data.frame(lowess(df$.overlap, df$.dissimilarity))
        p <- p + geom_line(aes(x = .data$x, y = .data$y), data = curve,
            colour = "#2c6fbb", linewidth = 1)
    }
    p + labs(x = overlap, y = dissimilarity)
}
