#' Plot transitions between states
#'
#' \code{plotTransitionMatrix()} counts the transitions between the states
#' of consecutive samples of each subject and draws them as a state by
#' state matrix.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param state.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains the state of each sample, such as a
#'   community type or the dominant taxon.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param as.relative \code{Logical scalar}. Shows the share of transitions
#'   from each state (\code{TRUE}) or the counts (\code{FALSE}).
#'   (Default: \code{TRUE})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the transitions into panels. Each
#'   transition takes the value of its earlier sample.
#'   (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Skips the transitions between
#'   consecutive samples of a subject that are further apart than this, in
#'   \code{time.col} units. For visits, \code{NULL} skips transitions
#'   across visits without a sample; otherwise \code{NULL} counts all of
#'   them. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_tile]{geom_tile}}.
#'
#' @details
#' Rows are the state before a transition and columns the state after it;
#' the diagonal holds the samples that stayed in their state. Samples with a
#' missing state break the chain: no transition is counted to or from them.
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
NULL

#' @rdname plotTransitionMatrix
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotTransitionMatrix", signature = c(x = "SummarizedExperiment"),
    function(x, state.col, time.col, group, as.relative = TRUE,
            facet.by = NULL, max.gap = NULL, ...) {
        .require_column(x, group)
        type <- .check_longitudinal(x, time.col, group)
        .require_column(x, state.col)
        .check_column(x, facet.by)
        if (!(isTRUE(as.relative) || isFALSE(as.relative))) {
            stop("'as.relative' must be TRUE or FALSE.", call. = FALSE)
        }
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
        if (is.null(max.gap) && type == "visit") {
            max.gap <- 1
        }
        df <- .break_gaps(df, time.col, group, max.gap)
        df <- .count_transitions(df, state.col, group, facet.by, as.relative)
        .plot_transition_matrix(df, state.col, facet.by, as.relative, ...)
    }
)

# Count transitions between consecutive samples of each subject, within a
# gap segment from .break_gaps(). Returns one row per from-to pair (and
# facet.by level), all pairs included, with columns from, to, n and value
# (count or share of the from row).
.count_transitions <- function(df, state.col, group, facet.by,
        as.relative) {
    state <- df[[state.col]]
    lv <- if (is.factor(state)) levels(state) else sort(unique(
        as.character(state[!is.na(state)])))
    i <- .next_in_subject(df[[group]])
    i <- i[df$.segment[i] == df$.segment[i + 1L] & !is.na(state[i]) &
        !is.na(state[i + 1L])]
    if (length(i) == 0L) {
        stop("No subject has two consecutive samples with a state.",
            call. = FALSE)
    }
    by <- if (is.null(facet.by)) "" else as.character(df[[facet.by]][i])
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
