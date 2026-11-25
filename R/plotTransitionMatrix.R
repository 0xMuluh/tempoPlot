#' Plot transitions between states
#'
#' \code{plotTransitionMatrix()} counts the transitions between the states
#' of consecutive samples of each subject and draws them as a state by
#' state matrix or as a graph of states joined by arrows.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param state.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains the state of each sample, such as a
#'   community type or the dominant taxon.
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
#' @param as.relative \code{Logical scalar}. Shows the share of transitions
#'   from each state (\code{TRUE}) or the counts (\code{FALSE}).
#'   (Default: \code{TRUE})
#'
#' @param layout \code{Character scalar}. Draws a state by state matrix
#'   (\code{"matrix"}) or the states on a circle joined by arrows
#'   (\code{"graph"}). (Default: \code{"matrix"})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the transitions into panels. Each
#'   transition takes the value of its later sample.
#'   (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Skips the transitions between
#'   consecutive samples of a subject that are further apart than this, in
#'   \code{time.col} units. For visits, \code{NULL} skips transitions
#'   across visits without a sample; otherwise \code{NULL} counts all of
#'   them. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_tile]{geom_tile}}, or to
#'   \code{\link[ggplot2:geom_segment]{geom_curve}} for the arrows of
#'   \code{layout = "graph"}.
#'
#' @details
#' Rows are the state before a transition and columns the state after it;
#' the diagonal holds the samples that stayed in their state. Samples with a
#' missing state break the chain: no transition is counted to or from them.
#'
#' In the graph, as in the state diagrams of Gajer et al. (2012,
#' doi:10.1126/scitranslmed.3003605, Fig. 5), the width of an arrow is the
#' share (or count) of transitions from one state to another. The size of a
#' state is the number of transitions from it, and its fill the share that
#' stayed in it.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotTransitionMatrix
#'
#' @examples
#' # Dominant genus of each sample: the four most common and "Other"
#' data("hitchip1006", package = "miaTime")
#' tse <- mia::addDominant(hitchip1006, name = "dominant", n = 4)
#' plotTransitionMatrix(tse, state.col = "dominant", time.col = "time",
#'     group = "subject")
#'
#' # Counts by sex
#' plotTransitionMatrix(tse, state.col = "dominant", time.col = "time",
#'     group = "subject", as.relative = FALSE, facet.by = "sex")
#'
#' # The same transitions as a graph
#' plotTransitionMatrix(tse, state.col = "dominant", time.col = "time",
#'     group = "subject", layout = "graph")
NULL

#' @rdname plotTransitionMatrix
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotTransitionMatrix", signature = c(x = "SummarizedExperiment"),
    function(x, state.col, time.col, group, as.relative = TRUE,
            layout = c("matrix", "graph"), facet.by = NULL, max.gap = NULL,
            ...) {
        layout <- match.arg(layout)
        .require_column(x, group)
        .check_longitudinal(x, time.col, group)
        .require_column(x, state.col)
        .check_column(x, facet.by)
        .check_flag(as.relative)
        state <- colData(x)[[state.col]]
        if (!(is.factor(state) || is.character(state) ||
                is.logical(state))) {
            stop("'state.col' must be a factor, character or logical ",
                "column.", call. = FALSE)
        }

        cols <- unique(c(time.col, group, state.col, facet.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df <- .order_time(df, time.col, group)
        df <- .break_gaps(df, time.col, group, max.gap)
        df <- .count_transitions(df, state.col, group, facet.by, as.relative)
        if (layout == "graph") {
            .plot_transition_graph(df, state.col, facet.by, as.relative, ...)
        } else {
            .plot_transition_matrix(df, state.col, facet.by, as.relative,
                ...)
        }
    }
)

# Count transitions between consecutive samples of each subject, within a
# gap segment from .break_gaps(). Returns one row per from-to pair (and
# facet.by level), all pairs included, with columns from, to, n and value
# (count or share of the from row).
.count_transitions <- function(df, state.col, group, facet.by,
        as.relative) {
    state <- df[[state.col]]
    lv <- if (is.factor(state)) levels(state) else .levels_of(state)
    i <- .next_in_subject(df[[group]], segment = df$.segment)
    i <- i[!is.na(state[i]) & !is.na(state[i + 1L])]
    if (length(i) == 0L) {
        stop("No subject has two consecutive samples with a state.",
            call. = FALSE)
    }
    by <- if (is.null(facet.by)) "" else as.character(df[[facet.by]][i + 1L])
    tr <- data.frame(
        from = factor(as.character(state[i]), levels = lv),
        to = factor(as.character(state[i + 1L]), levels = lv),
        by = factor(by))
    tab <- as.data.frame(table(tr), responseName = "n")
    tab$value <- tab$n
    if (as.relative) {
        total <- ave(tab$n, tab$from, tab$by, FUN = sum)
        tab$value <- ifelse(total > 0, tab$n / total, NA_real_)
    }
    if (!is.null(facet.by)) {
        names(tab)[names(tab) == "by"] <- facet.by
    }
    tab
}

#' @importFrom ggplot2 ggplot aes .data geom_tile geom_text labs
#'   facet_wrap vars scale_y_discrete theme element_text
.plot_transition_matrix <- function(df, state.col, facet.by, as.relative,
        ...) {
    label <- if (as.relative) sprintf("%.2f", df$value) else df$n
    df$.label <- ifelse(df$n > 0, label, "")
    # White labels on the dark lower half of the sequential scale.
    light <- df$value >= mean(range(df$value, na.rm = TRUE))
    df$.text <- ifelse(light %in% TRUE, "black", "white")
    p <- ggplot(df, aes(x = .data$to, y = .data$from)) +
        geom_tile(aes(fill = .data$value), ...) +
        geom_text(aes(label = .data$.label, colour = I(.data$.text)),
            size = 3) +
        scale_y_discrete(limits = rev) +
        .scale_sequential("fill") +
        labs(x = paste(state.col, "after"), y = paste(state.col, "before"),
            fill = if (as.relative) "share" else "count") +
        .theme_tempoPlot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    p
}

# The states of the transition table df on a circle, clockwise from the
# top: an arrow for each transition between two states, its width the
# value; each state sized by the transitions from it and filled by the
# share that stayed.
#' @importFrom ggplot2 ggplot aes .data geom_curve geom_point geom_text
#'   labs facet_wrap vars theme element_blank arrow unit coord_equal
#'   scale_linewidth_continuous scale_size_area scale_x_continuous
#'   scale_y_continuous expansion
#' @importFrom utils modifyList
.plot_transition_graph <- function(df, state.col, facet.by, as.relative,
        ...) {
    lv <- levels(df$from)
    angle <- pi / 2 - 2 * pi * (seq_along(lv) - 1L) / length(lv)
    by <- if (is.null(facet.by)) "" else df[[facet.by]]
    df$.out <- ave(df$n, df$from, by, FUN = sum)
    i <- as.integer(df$from)
    j <- as.integer(df$to)
    df$.x <- cos(angle[i])
    df$.y <- sin(angle[i])
    nodes <- df[i == j & df$.out > 0, , drop = FALSE]
    nodes$.stay <- nodes$n / nodes$.out
    edges <- df[i != j & df$n > 0, , drop = FALSE]
    # Arrows start and end short of the states they join.
    dx <- cos(angle[j[i != j & df$n > 0]]) - edges$.x
    dy <- sin(angle[j[i != j & df$n > 0]]) - edges$.y
    gap <- 0.15 / sqrt(dx^2 + dy^2)
    edges$.xend <- edges$.x + dx * (1 - gap)
    edges$.yend <- edges$.y + dy * (1 - gap)
    edges$.x <- edges$.x + dx * gap
    edges$.y <- edges$.y + dy * gap
    # Names outside the circle, aligned away from it.
    side <- round(cos(angle), 1L)
    labels <- data.frame(.label = vapply(strwrap(lv, 18L, simplify = FALSE),
            paste, character(1L), collapse = "\n"), .x = 1.2 * cos(angle),
        .y = 1.2 * sin(angle), .h = (1 - sign(side)) / 2,
        .v = ifelse(side == 0, (1 - sign(sin(angle))) / 2, 0.5))
    p <- ggplot() +
        do.call(geom_curve, modifyList(list(aes(x = .data$.x,
                y = .data$.y, xend = .data$.xend, yend = .data$.yend,
                linewidth = .data$value), data = edges, curvature = 0.2,
            colour = "grey40", arrow = arrow(length = unit(0.2, "cm"),
                type = "closed")), list(...))) +
        geom_point(aes(x = .data$.x, y = .data$.y, size = .data$.out,
            fill = .data$.stay), data = nodes, shape = 21,
            colour = "grey20") +
        geom_text(aes(x = .data$.x, y = .data$.y, label = .data$.label,
            hjust = .data$.h, vjust = .data$.v), data = labels, size = 3) +
        scale_linewidth_continuous(range = c(0.3, 2)) +
        scale_size_area(max_size = 10) +
        .scale_sequential("fill", limits = c(0, 1)) +
        scale_x_continuous(expand = expansion(mult = 0.6)) +
        scale_y_continuous(expand = expansion(mult = 0.15)) +
        coord_equal() +
        labs(linewidth = if (as.relative) "share" else "count",
            size = "transitions from", fill = "stayed",
            title = state.col) +
        .theme_tempoPlot() +
        theme(axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), panel.grid = element_blank(),
            legend.box = "vertical")
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    p
}
