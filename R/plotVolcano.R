#' Plot effects against their significance
#'
#' \code{plotVolcano()} draws the effect of each feature against
#' \eqn{-\log_{10}} of its p-value and labels the strongest hits.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object, with the results in \code{rowData(x)}, or a \code{data.frame}
#'   of results with one row per feature.
#'
#' @param effect \code{Character scalar}. Specifies the column holding the
#'   effects, such as log fold changes or slopes.
#'
#' @param p.col \code{Character scalar}. Specifies the column holding the
#'   p-values or adjusted p-values.
#'
#' @param threshold \code{Numeric scalar}. Features with \code{p.col} below
#'   this are significant. (Default: \code{0.05})
#'
#' @param label.top \code{Integer scalar}. Labels this many significant
#'   features with the smallest p-values. (Default: \code{10})
#'
#' @param feature.col \code{Character scalar}. Specifies the column of
#'   feature names that label the points, of \code{rowData(x)} or of the
#'   \code{data.frame}. If \code{NULL}, uses the row names.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Features with a missing effect or p-value are dropped. The dashed line
#' marks \code{threshold}.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotVolcano
#'
#' @examples
#' # Differences between men and women at the first sample of each subject
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006[, hitchip1006$time == 0 & !is.na(hitchip1006$sex)]
#' tse <- mia::transformAssay(tse, method = "clr", pseudocount = 1)
#' mat <- SummarizedExperiment::assay(tse, "clr")
#' male <- tse$sex == "male"
#' res <- data.frame(
#'     genus = rownames(tse),
#'     difference = rowMeans(mat[, male]) - rowMeans(mat[, !male]),
#'     p = apply(mat, 1, function(v) wilcox.test(v[male], v[!male])$p.value))
#' res$q <- p.adjust(res$p, method = "fdr")
#' plotVolcano(res, effect = "difference", p.col = "q", label.top = 5,
#'     feature.col = "genus")
NULL

#' @rdname plotVolcano
#' @importFrom SummarizedExperiment rowData
#' @importFrom methods setMethod
#' @export
setMethod("plotVolcano", signature = c(x = "SummarizedExperiment"),
    function(x, effect, p.col, threshold = 0.05, label.top = 10,
            feature.col = NULL, ...) {
        .check_column(x, effect, where = "rowData")
        .check_column(x, p.col, where = "rowData")
        .check_column(x, feature.col, where = "rowData")
        cols <- unique(c(effect, p.col, feature.col))
        df <- as.data.frame(rowData(x)[, cols, drop = FALSE],
            optional = TRUE)
        plotVolcano(df, effect, p.col, threshold, label.top, feature.col,
            ...)
    }
)

#' @rdname plotVolcano
#' @importFrom methods setMethod
#' @export
setMethod("plotVolcano", signature = c(x = "data.frame"),
    function(x, effect, p.col, threshold = 0.05, label.top = 10,
            feature.col = NULL, ...) {
        for (arg in c("effect", "p.col", "feature.col")) {
            col <- get(arg)
            if (arg == "feature.col" && is.null(col)) next
            if (!.is_one_name(col) || !col %in% names(x)) {
                stop("'", arg, "' must be one column name of 'x'.",
                    call. = FALSE)
            }
            if (arg != "feature.col" && !is.numeric(x[[col]])) {
                stop("'", arg, "' must be a numeric column.", call. = FALSE)
            }
        }
        .check_number(threshold, 0, 1, open = TRUE)
        .check_number(label.top, 0)
        p <- x[[p.col]]
        if (any(p < 0 | p > 1, na.rm = TRUE)) {
            stop("'p.col' must contain values between 0 and 1.",
                call. = FALSE)
        }

        df <- data.frame(
            .feature = if (is.null(feature.col)) rownames(x) else
                as.character(x[[feature.col]]),
            effect = x[[effect]], p = p)
        df <- df[!is.na(df$effect) & !is.na(df$p), , drop = FALSE]
        if (nrow(df) == 0L) {
            stop("No feature has both an effect and a p-value.",
                call. = FALSE)
        }
        df$significant <- df$p < threshold
        hits <- which(df$significant)
        hits <- hits[order(df$p[hits])][seq_len(min(label.top,
            length(hits)))]
        df$.label <- ""
        df$.label[hits] <- df$.feature[hits]
        .plot_volcano(df, effect, p.col, threshold, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_point geom_hline geom_vline labs
#' @importFrom ggrepel geom_text_repel
.plot_volcano <- function(df, effect, p.col, threshold, ...) {
    ggplot(df, aes(x = .data$effect, y = -log10(.data$p))) +
        geom_vline(xintercept = 0, colour = "grey60") +
        geom_hline(yintercept = -log10(threshold), linetype = "dashed",
            colour = "grey40") +
        geom_point(aes(colour = .data$significant), ...) +
        geom_text_repel(aes(label = .data$.label), size = 3,
            max.overlaps = Inf, min.segment.length = 0) +
        .scale_highlight("colour") +
        labs(x = effect, y = paste0("-log10(", p.col, ")"),
            colour = paste(p.col, "<", threshold)) +
        .theme_tempoPlot()
}
