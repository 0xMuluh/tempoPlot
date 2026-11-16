#' Plot a network of features or samples
#'
#' \code{plotNetwork()} draws a node-edge graph from pairs stored in
#' \code{x}: feature pairs (\code{rowPair}), sample pairs (\code{colPair})
#' or an edge table in \code{metadata(x)}, such as transitions between
#' states or associations estimated in time windows.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object; for \code{rowPair} and \code{colPair} a
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#'
#' @param edges \code{Character scalar}. Specifies the edges: the name of a
#'   \code{rowPair} or \code{colPair} of \code{x}, or of a \code{data.frame}
#'   in \code{metadata(x)} with the columns \code{from} and \code{to},
#'   optionally \code{weight} and the \code{split.by} column.
#'
#' @param split.by \code{Character scalar}. Specifies a column of the edge
#'   table that splits the edges into panels, for example time windows,
#'   drawn with the same node positions. (Default: \code{NULL})
#'
#' @param diff \code{Logical scalar}. With a \code{split.by} of two levels,
#'   draws one network that marks the edges only in the first level
#'   (lost), only in the second (gained) or in both (kept).
#'   (Default: \code{FALSE})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_segment]{geom_segment}}.
#'
#' @details
#' Requires igraph for the layout (Kamada-Kawai, which gives the same
#' positions on every run). Edges are coloured by the sign of
#' \code{weight} and their width follows its absolute value. Node names
#' are drawn for up to 40 nodes.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotNetwork
#'
#' @examples
#' # Correlations above 0.4 between 15 genera, in women and men
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006[, hitchip1006$time == 0 & !is.na(hitchip1006$sex)]
#' tse <- mia::transformAssay(tse, method = "clr", pseudocount = 1)
#' top <- names(sort(rowMeans(SummarizedExperiment::assay(tse, "counts")),
#'     decreasing = TRUE))[1:15]
#' edges <- do.call(rbind, lapply(c("female", "male"), function(s) {
#'     r <- cor(t(SummarizedExperiment::assay(tse, "clr")[top,
#'         tse$sex == s]))
#'     i <- which(abs(r) > 0.4 & upper.tri(r), arr.ind = TRUE)
#'     data.frame(from = top[i[, 1]], to = top[i[, 2]], weight = r[i],
#'         sex = s)
#' }))
#' S4Vectors::metadata(tse)$correlations <- edges
#' plotNetwork(tse, edges = "correlations", split.by = "sex")
#'
#' # Edges gained and lost from women to men
#' plotNetwork(tse, edges = "correlations", split.by = "sex", diff = TRUE)
NULL

#' @rdname plotNetwork
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotNetwork", signature = c(x = "SummarizedExperiment"),
    function(x, edges, split.by = NULL, diff = FALSE, ...) {
        if (!requireNamespace("igraph", quietly = TRUE)) {
            stop("plotNetwork() requires the igraph package.",
                call. = FALSE)
        }
        if (!(isTRUE(diff) || isFALSE(diff))) {
            stop("'diff' must be TRUE or FALSE.", call. = FALSE)
        }
        df <- .edge_table(x, edges)
        if (!is.null(split.by)) {
            if (!.is_one_name(split.by) || !split.by %in% names(df)) {
                stop("'split.by' must be one column of the edge table.",
                    call. = FALSE)
            }
            df$.split <- factor(df[[split.by]])
        }
        if (diff) {
            if (is.null(split.by) || nlevels(df$.split) != 2L) {
                stop("'diff' needs a 'split.by' with two levels.",
                    call. = FALSE)
            }
            df <- .edge_diff(df)
        }
        nodes <- unique(c(df$from, df$to))
        g <- igraph::graph_from_data_frame(
            unique(df[c("from", "to")]), vertices = nodes)
        pos <- igraph::layout_with_kk(g)
        pos <- data.frame(.node = nodes, .x = pos[, 1L], .y = pos[, 2L])
        df <- cbind(df,
            .x = pos$.x[match(df$from, nodes)],
            .y = pos$.y[match(df$from, nodes)],
            .xend = pos$.x[match(df$to, nodes)],
            .yend = pos$.y[match(df$to, nodes)])
        .plot_network(df, pos, split.by, diff, ...)
    }
)

# Edge table with from, to and weight (1 if absent), as characters.
#' @importFrom S4Vectors metadata
.edge_table <- function(x, edges) {
    if (!.is_one_name(edges)) {
        stop("'edges' must be a single name.", call. = FALSE)
    }
    if (inherits(x, "SingleCellExperiment")) {
        for (kind in c("row", "col")) {
            nms <- if (kind == "row")
                SingleCellExperiment::rowPairNames(x) else
                SingleCellExperiment::colPairNames(x)
            if (edges %in% nms) {
                h <- if (kind == "row")
                    SingleCellExperiment::rowPair(x, edges) else
                    SingleCellExperiment::colPair(x, edges)
                ids <- if (kind == "row") rownames(x) else colnames(x)
                ids <- ids %||% as.character(seq_along(ids))
                w <- S4Vectors::mcols(h)$value
                return(data.frame(from = ids[S4Vectors::from(h)],
                    to = ids[S4Vectors::to(h)],
                    weight = if (is.null(w)) 1 else as.numeric(w)))
            }
        }
    }
    df <- metadata(x)[[edges]]
    if (!is.data.frame(df) || !all(c("from", "to") %in% names(df))) {
        stop("'edges' must name a rowPair, a colPair or a data.frame ",
            "with the columns from and to in metadata(x).", call. = FALSE)
    }
    df <- as.data.frame(df)
    df$from <- as.character(df$from)
    df$to <- as.character(df$to)
    if (!"weight" %in% names(df)) {
        df$weight <- 1
    }
    df
}

# One table of the edges of two levels: lost, gained or kept.
.edge_diff <- function(df) {
    lv <- levels(df$.split)
    key <- paste(df$from, df$to, sep = "\r")
    in1 <- key %in% key[df$.split == lv[1L]]
    in2 <- key %in% key[df$.split == lv[2L]]
    df$.change <- ifelse(in1 & in2, "kept", ifelse(in1, paste("only",
        lv[1L]), paste("only", lv[2L])))
    df$.change <- factor(df$.change, levels = c(paste("only", lv[1L]),
        paste("only", lv[2L]), "kept"))
    df <- df[!duplicated(key) | df$.change != "kept", , drop = FALSE]
    df[!duplicated(paste(df$from, df$to, df$.change)), , drop = FALSE]
}

#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes .data geom_segment geom_point
#'   labs facet_wrap vars scale_colour_manual theme_void
#'   scale_linewidth_continuous scale_x_continuous scale_y_continuous
#'   expansion
#' @importFrom stats setNames
.plot_network <- function(df, pos, split.by, diff, ...) {
    p <- ggplot()
    if (diff) {
        lv <- levels(df$.change)
        p <- p + geom_segment(aes(x = .data$.x, y = .data$.y,
                xend = .data$.xend, yend = .data$.yend,
                colour = .data$.change), data = df, ...) +
            scale_colour_manual(values = setNames(
                c("#b2182b", "#2166ac", "grey60"), lv)) +
            labs(colour = NULL)
    } else {
        df$.sign <- factor(ifelse(df$weight < 0, "negative", "positive"),
            levels = c("positive", "negative"))
        p <- p + geom_segment(aes(x = .data$.x, y = .data$.y,
                xend = .data$.xend, yend = .data$.yend,
                colour = .data$.sign, linewidth = abs(.data$weight)),
                data = df, ...) +
            scale_colour_manual(values = c(positive = "#2166ac",
                negative = "#b2182b"), drop = FALSE) +
            scale_linewidth_continuous(range = c(0.3, 1.5)) +
            labs(colour = "weight", linewidth = "|weight|")
        if (!is.null(split.by)) {
            p <- p + facet_wrap(vars(.data$.split))
        }
    }
    p <- p + geom_point(aes(x = .data$.x, y = .data$.y), data = pos,
        size = 2.5, colour = "grey20")
    if (nrow(pos) <= 40L) {
        p <- p + geom_text_repel(aes(x = .data$.x, y = .data$.y,
            label = .data$.node), data = pos, size = 2.8,
            max.overlaps = Inf)
    }
    # Room for node names at the edges of the panels.
    p + scale_x_continuous(expand = expansion(mult = 0.2)) +
        scale_y_continuous(expand = expansion(mult = 0.1)) +
        theme_void()
}
