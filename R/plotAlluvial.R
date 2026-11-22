#' Plot flows of subjects between states over time
#'
#' \code{plotAlluvial()} draws the state of each subject at each time
#' point as stacked strata and the subjects moving between states as flows.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param state.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains the state of each sample, such as a
#'   community type or the dominant taxon.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains the time points shared by the subjects:
#'   numeric, \code{Date}, \code{POSIXct} or, best, a factor of visits.
#'   Each distinct time is one axis.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param min.flow \code{Integer scalar}. Hides the paths through all time
#'   points followed by fewer subjects than this. (Default: \code{1})
#'
#' @param missing \code{Character scalar}. Draws time points without a
#'   sample or a state as a \code{"missing"} stratum (\code{"show"}) or
#'   removes the subjects that have them (\code{"drop"}).
#'   (Default: \code{"show"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggalluvial:geom_flow]{ggalluvial::geom_flow}}.
#'
#' @details
#' Requires the ggalluvial package. Subjects that follow the same states
#' through all time points form one flow.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotAlluvial
#'
#' @examples
#' # Dominant genus at the first three samples of each subject
#' data("hitchip1006", package = "miaTime")
#' tse <- mia::addDominant(hitchip1006, name = "dominant", n = 4)
#' tse$visit <- ave(tse$time, tse$subject, FUN = rank)
#' tse <- tse[, tse$visit <= 3 &
#'     tse$subject %in% names(which(table(tse$subject) >= 2))]
#' tse$visit <- factor(tse$visit)
#' plotAlluvial(tse, state.col = "dominant", time.col = "visit",
#'     group = "subject")
NULL

#' @rdname plotAlluvial
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotAlluvial", signature = c(x = "SummarizedExperiment"),
    function(x, state.col, time.col, group, min.flow = 1,
            missing = c("show", "drop"), ...) {
        missing <- match.arg(missing)
        if (!requireNamespace("ggalluvial", quietly = TRUE)) {
            stop("plotAlluvial() requires the ggalluvial package.",
                call. = FALSE)
        }
        .require_column(x, group)
        .check_longitudinal(x, time.col, group)
        .require_column(x, state.col)
        .check_number(min.flow, 1)

        time <- colData(x)[[time.col]]
        axes <- if (is.factor(time)) levels(droplevels(time)) else
            as.character(sort(unique(time)))
        state <- colData(x)[[state.col]]
        states <- if (is.factor(state)) levels(droplevels(state)) else
            sort(unique(as.character(state[!is.na(state)])))
        subject <- as.character(colData(x)[[group]])
        wide <- matrix(NA_character_, length(unique(subject)), length(axes),
            dimnames = list(unique(subject), axes))
        wide[cbind(subject, as.character(time))] <- as.character(state)
        df <- .alluvial_flows(wide, states, missing, min.flow)
        .plot_alluvial(df, state.col, time.col, ...)
    }
)

# One row per distinct path of states through the time points (columns of
# wide, one row per subject) and time point, weighted by its number of
# subjects. Gaps become the state "missing" or drop their subject.
.alluvial_flows <- function(wide, states, missing, min.flow) {
    if (missing == "drop") {
        gapped <- rowSums(is.na(wide)) > 0
        if (all(gapped)) {
            stop("All subjects have missing time points.", call. = FALSE)
        }
        if (any(gapped)) {
            message("Removed ", sum(gapped),
                " subjects with missing time points.")
        }
        wide <- wide[!gapped, , drop = FALSE]
    } else if (anyNA(wide)) {
        wide[is.na(wide)] <- "missing"
        states <- c(setdiff(states, "missing"), "missing")
    }
    path <- apply(wide, 1L, paste, collapse = "\r")
    n <- table(path)
    keep <- names(n)[n >= min.flow]
    if (length(keep) == 0L) {
        stop("No path is followed by 'min.flow' subjects.", call. = FALSE)
    }
    if (length(keep) < length(n)) {
        message("Hid ", sum(n[n < min.flow]), " subjects on paths ",
            "followed by fewer than ", min.flow, " subjects.")
    }
    first <- wide[match(keep, path), , drop = FALSE]
    axes <- colnames(wide)
    data.frame(
        .flow = rep(seq_along(keep), times = length(axes)),
        .time = factor(rep(axes, each = length(keep)), levels = axes),
        .state = factor(as.vector(first), levels = states),
        .n = rep(as.vector(n[keep]), times = length(axes)))
}

#' @importFrom ggplot2 ggplot aes .data labs
.plot_alluvial <- function(df, state.col, time.col, ...) {
    ggplot(df, aes(x = .data$.time, stratum = .data$.state,
            alluvium = .data$.flow, y = .data$.n, fill = .data$.state)) +
        ggalluvial::geom_flow(...) +
        ggalluvial::geom_stratum(colour = "white", linewidth = 0.3) +
        .scale_categorical("fill", grey = "missing") +
        labs(x = time.col, y = "subjects", fill = state.col) +
        .theme_tempoPlot()
}
