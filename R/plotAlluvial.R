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
#'   \code{colData(x)} that contains the time points shared by the subjects,
#'   best a factor of visits. Each distinct time is one axis.
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
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        .time_type(x, time.col)
        .check_column(x, state.col)
        if (is.null(state.col)) {
            stop("'state.col' must be one column name of colData(x).",
                call. = FALSE)
        }
        if (!(is.numeric(min.flow) && length(min.flow) == 1L &&
                min.flow >= 1)) {
            stop("'min.flow' must be a single number of at least 1.",
                call. = FALSE)
        }

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

        if (missing == "drop") {
            gapped <- rowSums(is.na(wide)) > 0
            if (all(gapped)) {
                stop("All subjects have missing time points.",
                    call. = FALSE)
            }
            if (any(gapped)) {
                message("Removed ", sum(gapped),
                    " subjects with missing time points.")
            }
            wide <- wide[!gapped, , drop = FALSE]
        } else {
            wide[is.na(wide)] <- "missing"
            if (any(wide == "missing")) {
                states <- c(setdiff(states, "missing"), "missing")
            }
        }

        # One flow per distinct path, weighted by its number of subjects.
        path <- apply(wide, 1L, paste, collapse = "\r")
        n <- table(path)
        keep <- names(n)[n >= min.flow]
        if (length(keep) == 0L) {
            stop("No path is followed by 'min.flow' subjects.",
                call. = FALSE)
        }
        if (length(keep) < length(n)) {
            message("Hid ", sum(n[n < min.flow]), " subjects on paths ",
                "followed by fewer than ", min.flow, " subjects.")
        }
        first <- wide[match(keep, path), , drop = FALSE]
        df <- data.frame(
            .flow = rep(seq_along(keep), times = length(axes)),
            .time = factor(rep(axes, each = length(keep)), levels = axes),
            .state = factor(as.vector(first), levels = states),
            .n = rep(as.vector(n[keep]), times = length(axes)))
        .plot_alluvial(df, state.col, time.col, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data labs scale_fill_manual
#' @importFrom grDevices hcl
.plot_alluvial <- function(df, state.col, time.col, ...) {
    states <- levels(df$.state)
    k <- length(setdiff(states, "missing"))
    colours <- hcl(h = seq(15, 375, length.out = k + 1L)[seq_len(k)],
        c = 100, l = 65)
    names(colours) <- setdiff(states, "missing")
    colours <- c(colours, missing = "grey80")
    ggplot(df, aes(x = .data$.time, stratum = .data$.state,
            alluvium = .data$.flow, y = .data$.n, fill = .data$.state)) +
        ggalluvial::geom_flow(...) +
        ggalluvial::geom_stratum(colour = "white", linewidth = 0.3) +
        scale_fill_manual(values = colours) +
        labs(x = time.col, y = "subjects", fill = state.col)
}
