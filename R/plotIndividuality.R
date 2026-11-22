#' Plot within- against between-subject dissimilarities
#'
#' \code{plotIndividuality()} compares the dissimilarities between samples
#' of the same subject with those between samples of different subjects.
#' The further apart the two distributions, the more individual the
#' communities.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param dissimilarity \code{Character scalar}. Specifies the
#'   sample by sample dissimilarity in \code{metadata(x)}, as stored by
#'   \code{\link[mia:getDissimilarity]{mia::addDissimilarity()}}.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param layout \code{Character scalar}. Draws density curves
#'   (\code{"density"}), box plots (\code{"box"}) or violins
#'   (\code{"violin"}). (Default: \code{"density"})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the pairs into panels. Only pairs of
#'   samples from the same group are kept. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_density]{geom_density}},
#'   \code{\link[ggplot2:geom_boxplot]{geom_boxplot}} or
#'   \code{\link[ggplot2:geom_violin]{geom_violin}}.
#'
#' @details
#' Each pair of samples is counted once. A subject with \eqn{k} samples
#' gives \eqn{k(k-1)/2} within-subject pairs, so subjects sampled often
#' weigh most in the within distribution.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotIndividuality
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 2))]
#' tse <- mia::addDissimilarity(tse, method = "bray", assay.type = "counts")
#' plotIndividuality(tse, dissimilarity = "bray", group = "subject")
#'
#' # By sex
#' plotIndividuality(tse, dissimilarity = "bray", group = "subject",
#'     layout = "box", facet.by = "sex")
NULL

#' @rdname plotIndividuality
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotIndividuality", signature = c(x = "SummarizedExperiment"),
    function(x, dissimilarity, group, layout = c("density", "box", "violin"),
            facet.by = NULL, ...) {
        layout <- match.arg(layout)
        .require_column(x, group)
        .check_column(x, facet.by)
        d <- .get_dissimilarity(x, dissimilarity)

        pairs <- which(lower.tri(d), arr.ind = TRUE)
        subject <- as.character(colData(x)[[group]])
        same <- subject[pairs[, 1L]] == subject[pairs[, 2L]]
        df <- data.frame(
            dissimilarity = d[pairs],
            pair = factor(ifelse(same, "within", "between"),
                levels = c("within", "between")))
        if (!is.null(facet.by)) {
            by <- as.character(colData(x)[[facet.by]])
            keep <- by[pairs[, 1L]] == by[pairs[, 2L]]
            keep <- !is.na(keep) & keep
            df[[facet.by]] <- by[pairs[, 1L]]
            df <- df[keep, , drop = FALSE]
        }
        df <- df[!is.na(df$dissimilarity), , drop = FALSE]
        if (!any(df$pair == "within")) {
            stop("No subject has more than one sample.", call. = FALSE)
        }
        .plot_individuality(df, dissimilarity, layout, facet.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_density geom_boxplot
#'   geom_violin labs facet_wrap vars
.plot_individuality <- function(df, dissimilarity, layout, facet.by, ...) {
    if (layout == "density") {
        p <- ggplot(df, aes(x = .data$dissimilarity,
                colour = .data$pair)) +
            geom_density(...) +
            .scale_categorical("colour") +
            labs(x = dissimilarity, colour = "pair")
    } else {
        geom <- if (layout == "box") geom_boxplot else geom_violin
        p <- ggplot(df, aes(x = .data$pair, y = .data$dissimilarity)) +
            geom(...) +
            labs(x = "pair", y = dissimilarity)
    }
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    p + .theme_tempoPlot()
}
