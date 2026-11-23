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
#' @param directed \code{Logical scalar}. Whether edges point from
#'   \code{from} to \code{to}, such as transitions between states. Directed
#'   edges are drawn as arrows, and \code{diff} tells a-b from b-a;
#'   otherwise they are the same edge. (Default: \code{FALSE})
#'
#' @param layout \code{Character scalar}. Draws nodes joined by edges
#'   (\code{"graph"}) or a matrix of tiles, \code{from} in rows and
#'   \code{to} in columns (\code{"matrix"}). (Default: \code{"graph"})
#'
#' @param p.col \code{Character scalar}. With \code{layout = "matrix"}, a
#'   numeric column of the edge table holding p-values or adjusted
#'   p-values; tiles below \code{threshold} are marked with a dot.
#'   (Default: \code{NULL})
#'
#' @param threshold \code{Numeric scalar}. Significance level for
#'   \code{p.col}. (Default: \code{0.05})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_segment]{geom_segment}}, or
#'   \code{\link[ggplot2:geom_tile]{geom_tile}} for the matrix.
#'
#' @details
#' Requires igraph for the layout (Kamada-Kawai, which gives the same
#' positions on every run). Edges are coloured by the sign of
#' \code{weight} and their width follows its absolute value. Self-loops
#' (\code{from} equal to \code{to}), such as staying in a state, are drawn
#' as rings around their node. With \code{diff}, an edge whose weight
#' changes sign is kept. Node names are drawn for up to 40 nodes.
#'
#' The matrix shows interaction strengths, such as those of a generalised
#' Lotka-Volterra model (Bucci et al. 2016, Fig. 4a,
#' \doi{10.1186/s13059-016-0980-6}). Tiles are filled by \code{weight},
#' blue below zero and red above, or by the change with \code{diff}.
#' Self-loops fill the diagonal; pairs without an edge stay blank. Nodes are
#' sorted by name in rows and columns, so a node keeps its place in every
#' panel. Undirected edges fill both a-b and b-a, so the matrix is
#' symmetric. The matrix does not need igraph.
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
#'
#' # All correlations as a matrix, with FDR q-values below 0.05 marked
#' pairs <- t(combn(top, 2))
#' edges <- do.call(rbind, lapply(c("female", "male"), function(s) {
#'     m <- SummarizedExperiment::assay(tse, "clr")[, tse$sex == s]
#'     test <- apply(pairs, 1, function(p) cor.test(m[p[1], ], m[p[2], ]))
#'     data.frame(from = pairs[, 1], to = pairs[, 2], sex = s,
#'         weight = vapply(test, `[[`, numeric(1), "estimate"),
#'         p = vapply(test, `[[`, numeric(1), "p.value"))
#' }))
#' edges$q <- p.adjust(edges$p, method = "fdr")
#' S4Vectors::metadata(tse)$all <- edges
#' plotNetwork(tse, edges = "all", split.by = "sex", layout = "matrix",
#'     p.col = "q")
#'
#' # Moves between the dominant genera of consecutive samples
#' tse <- mia::addDominant(hitchip1006, name = "dominant")
#' cd <- as.data.frame(SummarizedExperiment::colData(tse))
#' cd <- cd[order(cd$subject, cd$time), ]
#' n <- nrow(cd)
#' same <- cd$subject[-1] == cd$subject[-n]
#' moves <- data.frame(from = cd$dominant[-n][same],
#'     to = cd$dominant[-1][same], weight = 1)
#' moves <- aggregate(weight ~ from + to, data = moves, FUN = sum)
#' S4Vectors::metadata(tse)$moves <- moves
#' plotNetwork(tse, edges = "moves", directed = TRUE)
NULL

#' @rdname plotNetwork
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotNetwork", signature = c(x = "SummarizedExperiment"),
    function(x, edges, split.by = NULL, diff = FALSE, directed = FALSE,
            layout = c("graph", "matrix"), p.col = NULL, threshold = 0.05,
            ...) {
        layout <- match.arg(layout)
        .check_flag(diff)
        .check_flag(directed)
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
            df <- .edge_diff(df, directed)
        }
        if (layout == "matrix") {
            return(.plot_edge_matrix(df, split.by, diff, directed, p.col,
                threshold, ...))
        }
        if (!is.null(p.col)) {
            stop("'p.col' needs layout = \"matrix\".", call. = FALSE)
        }
        if (!requireNamespace("igraph", quietly = TRUE)) {
            stop("plotNetwork() requires the igraph package.",
                call. = FALSE)
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
        if (directed) {
            # Arrows stop short of the target node.
            df$.xend <- df$.x + 0.9 * (df$.xend - df$.x)
            df$.yend <- df$.y + 0.9 * (df$.yend - df$.y)
        }
        .plot_network(df, pos, split.by, diff, directed, ...)
    }
)

# Edge table with from, to and weight (1 if absent), as characters.
#' @importFrom S4Vectors metadata mcols from to nnode
#' @importFrom SingleCellExperiment rowPair colPair rowPairNames
#'   colPairNames
.edge_table <- function(x, edges) {
    if (!.is_one_name(edges)) {
        stop("'edges' must be a single name.", call. = FALSE)
    }
    if (inherits(x, "SingleCellExperiment") &&
            edges %in% c(rowPairNames(x), colPairNames(x))) {
        row <- edges %in% rowPairNames(x)
        h <- if (row) rowPair(x, edges) else colPair(x, edges)
        ids <- (if (row) rownames(x) else colnames(x)) %||%
            as.character(seq_len(nnode(h)))
        w <- mcols(h)$value
        return(data.frame(from = ids[from(h)], to = ids[to(h)],
            weight = if (is.null(w)) 1 else as.numeric(w)))
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

# One table of the edges of two levels: lost, gained or kept. Undirected
# edges a-b and b-a are the same edge.
.edge_diff <- function(df, directed) {
    lv <- levels(df$.split)
    key <- if (directed) paste(df$from, df$to, sep = "\r") else
        paste(pmin(df$from, df$to), pmax(df$from, df$to), sep = "\r")
    in1 <- key %in% key[df$.split == lv[1L]]
    in2 <- key %in% key[df$.split == lv[2L]]
    df$.change <- ifelse(in1 & in2, "kept", ifelse(in1, paste("only",
        lv[1L]), paste("only", lv[2L])))
    df$.change <- factor(df$.change, levels = c(paste("only", lv[1L]),
        paste("only", lv[2L]), "kept"))
    df[!duplicated(key), , drop = FALSE]
}

#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes .data geom_segment geom_point
#'   labs facet_wrap vars theme element_blank
#'   scale_linewidth_continuous scale_x_continuous scale_y_continuous
#'   expansion arrow unit
#' @importFrom utils modifyList
.plot_network <- function(df, pos, split.by, diff, directed, ...) {
    colour <- if (diff) ".change" else ".sign"
    df$.sign <- factor(ifelse(df$weight < 0, "negative", "positive"),
        levels = c("negative", "positive"))
    loop <- df$from == df$to
    edge <- aes(x = .data$.x, y = .data$.y, xend = .data$.xend,
        yend = .data$.yend, colour = .data[[colour]],
        linewidth = if (!diff) abs(.data$weight))
    p <- ggplot() +
        do.call(geom_segment, modifyList(list(edge,
            data = df[!loop, , drop = FALSE], arrow = if (directed)
                arrow(length = unit(0.2, "cm"), type = "closed")),
            list(...)))
    if (any(loop)) {
        # Self-loops as rings around their node.
        p <- p + geom_point(aes(x = .data$.x, y = .data$.y,
                colour = .data[[colour]]), data = df[loop, , drop = FALSE],
            shape = 1, size = 5, stroke = 1)
    }
    p <- p + .scale_categorical("colour", drop = FALSE) +
        scale_linewidth_continuous(range = c(0.3, 1.5)) +
        labs(colour = if (diff) split.by else "weight",
            linewidth = "|weight|")
    if (!diff && !is.null(split.by)) {
        p <- p + facet_wrap(vars(.data$.split))
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
        .theme_tempoPlot() +
        theme(axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), panel.grid = element_blank())
}

# Tiles of the edges, from in rows and to in columns, nodes sorted by name.
# Undirected edges fill both a-b and b-a; a table that already lists both
# orientations collapses to one tile each.
#' @importFrom ggplot2 ggplot aes .data geom_tile geom_point labs
#'   facet_wrap vars scale_x_discrete scale_y_discrete coord_equal theme
#'   element_text
.plot_edge_matrix <- function(df, split.by, diff, directed, p.col,
        threshold, ...) {
    if (!is.null(p.col) && !(.is_one_name(p.col) && p.col %in% names(df) &&
            is.numeric(df[[p.col]]))) {
        stop("'p.col' must be one numeric column of the edge table.",
            call. = FALSE)
    }
    .check_number(threshold, 0, 1, open = TRUE)
    if (!directed) {
        flip <- df[df$from != df$to, , drop = FALSE]
        flip[c("from", "to")] <- flip[c("to", "from")]
        df <- unique(rbind(df, flip))
    }
    if (anyDuplicated(df[intersect(c(".split", "from", "to"), names(df))])) {
        stop("'edges' must have one row per pair of nodes",
            if (!is.null(split.by)) " in each 'split.by' level", ".",
            call. = FALSE)
    }
    nodes <- sort(unique(c(df$from, df$to)))
    fill <- if (diff) ".change" else "weight"
    p <- ggplot(df, aes(x = .data$to, y = .data$from)) +
        geom_tile(aes(fill = .data[[fill]]), ...) +
        (if (diff) .scale_categorical("fill", drop = FALSE) else
            .scale_diverging("fill")) +
        scale_x_discrete(limits = nodes) +
        scale_y_discrete(limits = rev(nodes)) +
        coord_equal() +
        labs(x = "to", y = "from", fill = if (diff) split.by else "weight") +
        .theme_tempoPlot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (!is.null(p.col)) {
        p <- p + geom_point(data = df[which(df[[p.col]] < threshold), ,
                drop = FALSE], size = 1) +
            labs(caption = paste("Dots:", p.col, "<", threshold))
    }
    if (!diff && !is.null(split.by)) {
        p <- p + facet_wrap(vars(.data$.split))
    }
    p
}
