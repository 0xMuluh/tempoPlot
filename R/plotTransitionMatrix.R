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
#' @param normalize \code{Character scalar}. Shows the share of transitions
#'   from each state (\code{"row"}) or the counts (\code{"none"}).
#'   (Default: \code{"row"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the transitions into panels. Each
#'   transition takes the value of its earlier sample.
#'   (Default: \code{NULL})
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
#'     group = "subject", normalize = "none", group.by = "sex")
NULL

#' @rdname plotTransitionMatrix
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotTransitionMatrix", signature = c(x = "SummarizedExperiment"),
    function(x, state.col, time.col, group, normalize = c("row", "none"),
            group.by = NULL, ...) {
        normalize <- match.arg(normalize)
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        .time_type(x, time.col)
        .check_column(x, state.col)
        .check_column(x, group.by)
        if (is.null(state.col)) {
            stop("'state.col' must be one column name of colData(x).",
                call. = FALSE)
        }
        state <- colData(x)[[state.col]]
        if (!(is.factor(state) || is.character(state) ||
                is.logical(state))) {
            stop("'state.col' must be a factor, character or logical ",
                "column.", call. = FALSE)
        }

        cols <- unique(c(time.col, group, state.col, group.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df <- .order_time(df, time.col, group)
        df <- .count_transitions(df, state.col, group, group.by, normalize)
        .plot_transition_matrix(df, state.col, group.by, normalize, ...)
    }
)

# Count transitions between consecutive samples of each subject. Returns
# one row per from-to pair (and group.by level), all pairs included, with
# columns from, to, n and value (count or share of the from row).
.count_transitions <- function(df, state.col, group, group.by, normalize) {
    state <- df[[state.col]]
    lv <- if (is.factor(state)) levels(state) else sort(unique(
        as.character(state[!is.na(state)])))
    n <- nrow(df)
    keep <- df[[group]][-1L] == df[[group]][-n] &
        !is.na(state[-1L]) & !is.na(state[-n])
    if (!any(keep)) {
        stop("No subject has two consecutive samples with a state.",
            call. = FALSE)
    }
    by <- if (is.null(group.by)) rep("", n - 1L) else
        as.character(df[[group.by]][-n])
    tr <- data.frame(
        from = factor(as.character(state[-n]), levels = lv)[keep],
        to = factor(as.character(state[-1L]), levels = lv)[keep],
        by = factor(by[keep]))
    tab <- as.data.frame(table(tr), responseName = "n")
    tab$value <- tab$n
    if (normalize == "row") {
        total <- ave(tab$n, tab$from, tab$by, FUN = sum)
        tab$value <- ifelse(total > 0, tab$n / total, NA_real_)
    }
    if (!is.null(group.by)) {
        names(tab)[names(tab) == "by"] <- group.by
    }
    tab
}

#' @importFrom ggplot2 ggplot aes .data geom_tile geom_text labs
#'   facet_wrap vars scale_y_discrete scale_fill_gradient theme
#'   element_text
.plot_transition_matrix <- function(df, state.col, group.by, normalize,
        ...) {
    label <- if (normalize == "row") round(df$value, 2) else df$n
    df$.label <- ifelse(df$n > 0, label, "")
    p <- ggplot(df, aes(x = .data$to, y = .data$from)) +
        geom_tile(aes(fill = .data$value), ...) +
        geom_text(aes(label = .data$.label), size = 3) +
        scale_y_discrete(limits = rev) +
        scale_fill_gradient(low = "grey95", high = "#5b9bd5",
            na.value = "white") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = paste(state.col, "after"), y = paste(state.col, "before"),
            fill = if (normalize == "row") "share" else "count")
    if (!is.null(group.by)) {
        p <- p + facet_wrap(vars(.data[[group.by]]))
    }
    p
}
