#' Plot time to an event
#'
#' \code{plotSurvival()} draws Kaplan-Meier curves of the time from each
#' subject's first sample to an event: dropout, the arrival or loss of a
#' feature, leaving the first state, or a clinical outcome.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param event \code{Character scalar}. The event: the last sample of a
#'   subject is before the latest sampling time in \code{x}
#'   (\code{"dropout"}); \code{feature} is present after being absent at
#'   the first sample (\code{"arrival"}) or absent after being present
#'   (\code{"loss"}); \code{state.col} differs from the first sample
#'   (\code{"state.exit"}); or \code{event.col} is \code{TRUE}
#'   (\code{"outcome"}). (Default: \code{"dropout"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the subjects into curves, taken from
#'   their first sample. (Default: \code{NULL})
#'
#' @param feature \code{Character scalar}. Specifies the feature for
#'   \code{event = "arrival"} and \code{"loss"}. A feature is present when
#'   its value in \code{assay.type} is above zero. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay for
#'   \code{feature}. (Default: \code{"counts"})
#'
#' @param state.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} with the state of each sample, for
#'   \code{event = "state.exit"}. (Default: \code{NULL})
#'
#' @param event.col \code{Character scalar}. Specifies the logical column
#'   from \code{colData(x)} that marks the samples where the outcome has
#'   occurred, for \code{event = "outcome"}. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_step]{geom_step}}.
#'
#' @details
#' Follow-up starts at the first sample of each subject. A subject without
#' the event is censored at its last sample, marked with a cross. For
#' \code{"arrival"} and \code{"loss"}, only subjects without and with the
#' feature at their first sample are at risk. Times are in \code{time.col}
#' units: days for \code{Date} times and visit steps for visits.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotSurvival
#'
#' @examples
#' # Time until the dominant genus changes, by sex
#' data("hitchip1006", package = "miaTime")
#' tse <- mia::addDominant(hitchip1006, name = "dominant")
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 2))]
#' plotSurvival(tse, time.col = "time", group = "subject",
#'     event = "state.exit", state.col = "dominant", group.by = "sex")
NULL

#' @rdname plotSurvival
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom methods setMethod
#' @export
setMethod("plotSurvival", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group,
            event = c("dropout", "arrival", "loss", "state.exit",
                "outcome"),
            group.by = NULL, feature = NULL, assay.type = "counts",
            state.col = NULL, event.col = NULL, ...) {
        event <- match.arg(event)
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        .time_type(x, time.col)
        .check_column(x, group.by)
        mark <- .event_marks(x, event, feature, assay.type, state.col,
            event.col)

        tnum <- .time_numeric(colData(x)[[time.col]])
        subject <- as.character(colData(x)[[group]])
        by <- if (is.null(group.by)) rep("", ncol(x)) else
            as.character(colData(x)[[group.by]])
        last.time <- max(tnum)
        rows <- lapply(split(seq_len(ncol(x)), subject), function(i) {
            i <- i[order(tnum[i])]
            t <- tnum[i] - tnum[i[1L]]
            if (event == "dropout") {
                hit <- if (tnum[i[length(i)]] < last.time) length(i)
            } else {
                m <- mark(i)
                if (is.na(m$at.risk) || !m$at.risk) {
                    return(NULL)
                }
                hit <- which(m$hit)[1L]
                if (is.na(hit)) hit <- NULL
            }
            data.frame(.by = by[i[1L]],
                .time = if (is.null(hit)) t[length(t)] else t[hit],
                .status = as.integer(!is.null(hit)))
        })
        df <- do.call(rbind, rows)
        if (is.null(df) || nrow(df) == 0L) {
            stop("No subject is at risk of the event.", call. = FALSE)
        }
        df <- .kaplan_meier(df)
        .plot_survival(df, time.col, group.by, event, ...)
    }
)

# Function of the ordered sample indices of one subject, returning whether
# the subject is at risk (at.risk) and at which samples the event has
# occurred (hit, the first sample never counts).
.event_marks <- function(x, event, feature, assay.type, state.col,
        event.col) {
    if (event == "dropout") {
        return(NULL)
    }
    if (event %in% c("arrival", "loss")) {
        if (!.is_one_name(assay.type) || !assay.type %in% assayNames(x)) {
            stop("'assay.type' must be one of assayNames(x).",
                call. = FALSE)
        }
        if (!.is_one_name(feature) || !feature %in% rownames(x)) {
            stop("'feature' must be one row name of 'x'.", call. = FALSE)
        }
        present <- assay(x, assay.type)[feature, ] > 0
        # The event state: present for arrival, absent for loss.
        target <- if (event == "arrival") present else !present
        return(function(i) list(at.risk = !target[i[1L]],
            hit = c(FALSE, target[i[-1L]])))
    }
    if (event == "state.exit") {
        .check_column(x, state.col)
        if (is.null(state.col)) {
            stop("'state.col' must be one column name of colData(x).",
                call. = FALSE)
        }
        state <- as.character(colData(x)[[state.col]])
        return(function(i) list(at.risk = !is.na(state[i[1L]]),
            hit = c(FALSE, !is.na(state[i[-1L]]) &
                state[i[-1L]] != state[i[1L]])))
    }
    .check_column(x, event.col)
    if (is.null(event.col) || !is.logical(colData(x)[[event.col]])) {
        stop("'event.col' must be one logical column of colData(x).",
            call. = FALSE)
    }
    outcome <- colData(x)[[event.col]]
    function(i) list(at.risk = !isTRUE(outcome[i[1L]]),
        hit = c(FALSE, outcome[i[-1L]] %in% TRUE))
}

# Kaplan-Meier estimate per .by: one row per step, with the share of
# subjects without the event (.surv) and censoring marks (.censored).
#' @importFrom survival survfit Surv
.kaplan_meier <- function(df) {
    res <- lapply(split(df, df$.by), function(d) {
        fit <- survfit(Surv(d$.time, d$.status) ~ 1)
        data.frame(.by = d$.by[1L],
            .time = c(0, fit$time), .surv = c(1, fit$surv),
            .censored = c(0L, fit$n.censor))
    })
    do.call(rbind, res)
}

#' @importFrom ggplot2 ggplot aes .data geom_step geom_point labs ylim
.plot_survival <- function(df, time.col, group.by, event, ...) {
    p <- ggplot(df, aes(x = .data$.time, y = .data$.surv,
            colour = if (!is.null(group.by)) .data$.by)) +
        geom_step(...) +
        geom_point(data = df[df$.censored > 0L, , drop = FALSE],
            shape = 3) +
        ylim(0, 1) +
        labs(x = paste(time.col, "from the first sample"),
            y = switch(event, dropout = "share still sampled",
                arrival = "share before arrival",
                loss = "share before loss",
                state.exit = "share in the first state",
                outcome = "share without the outcome"),
            colour = group.by)
    p
}
