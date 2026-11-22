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
#'   its value in \code{assay.type} is above \code{detection}.
#'   (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay for
#'   \code{feature}. (Default: \code{"counts"})
#'
#' @param detection \code{Numeric scalar}. Specifies the value above which
#'   \code{feature} is present. (Default: \code{0})
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
#' the event is censored at its last sample, marked with a cross. The band
#' is the 95\% confidence interval of the Kaplan-Meier estimate. Events
#' are only seen at sampling times, so the event time is the first sample
#' at which the event is observed, an upper bound of the true time.
#'
#' A subject drops out when its last sample precedes the last sampling
#' time of the study, the latest time in \code{x}. A subject that joins
#' late and is sampled until the end is not a dropout. For
#' \code{"arrival"} and \code{"loss"}, only subjects without and with the
#' feature at their first sample are at risk. Times are in \code{time.col}
#' units: days for \code{Date} times and visit steps for visits.
#'
#' @references
#' Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete
#' observations. \emph{Journal of the American Statistical Association}
#' 53(282):457-481.
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
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotSurvival", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group,
            event = c("dropout", "arrival", "loss", "state.exit",
                "outcome"),
            group.by = NULL, feature = NULL, assay.type = "counts",
            detection = 0, state.col = NULL, event.col = NULL, ...) {
        event <- match.arg(event)
        .require_column(x, group)
        .check_longitudinal(x, time.col, group)
        .check_column(x, group.by)
        mark <- .event_marks(x, event, feature, assay.type, detection,
            state.col, event.col)

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
.event_marks <- function(x, event, feature, assay.type, detection,
        state.col, event.col) {
    if (event == "dropout") {
        return(NULL)
    }
    if (event %in% c("arrival", "loss")) {
        if (!.is_one_name(feature) || !feature %in% rownames(x)) {
            stop("'feature' must be one row name of 'x'.", call. = FALSE)
        }
        present <- .present(.get_assay(x, assay.type)[feature, ], detection)
        # The event state: present for arrival, absent for loss.
        target <- if (event == "arrival") present else !present
        return(function(i) list(at.risk = !target[i[1L]],
            hit = c(FALSE, target[i[-1L]])))
    }
    if (event == "state.exit") {
        .require_column(x, state.col)
        state <- as.character(colData(x)[[state.col]])
        return(function(i) list(at.risk = !is.na(state[i[1L]]),
            hit = c(FALSE, !is.na(state[i[-1L]]) &
                state[i[-1L]] != state[i[1L]])))
    }
    .require_column(x, event.col)
    if (!is.logical(colData(x)[[event.col]])) {
        stop("'event.col' must be one logical column of colData(x).",
            call. = FALSE)
    }
    outcome <- colData(x)[[event.col]]
    function(i) list(at.risk = !isTRUE(outcome[i[1L]]),
        hit = c(FALSE, outcome[i[-1L]] %in% TRUE))
}

# Kaplan-Meier estimate per .by: one row per step, with the share of
# subjects without the event (.surv), its confidence interval (.lower,
# .upper) until the next step (.next) and censoring marks (.censored).
#' @importFrom survival survfit Surv
.kaplan_meier <- function(df) {
    res <- lapply(split(df, df$.by), function(d) {
        fit <- survfit(Surv(d$.time, d$.status) ~ 1)
        time <- c(0, fit$time)
        data.frame(.by = d$.by[1L], .time = time,
            .next = c(fit$time, time[length(time)]),
            .surv = c(1, fit$surv), .lower = c(1, fit$lower),
            .upper = c(1, fit$upper), .censored = c(0L, fit$n.censor))
    })
    do.call(rbind, res)
}

#' @importFrom ggplot2 ggplot aes .data geom_step geom_point geom_rect labs
#'   ylim
.plot_survival <- function(df, time.col, group.by, event, ...) {
    by <- !is.null(group.by)
    ggplot(df, aes(x = .data$.time, y = .data$.surv,
            colour = if (by) .data$.by)) +
        geom_rect(aes(xmin = .data$.time, xmax = .data$.next,
            ymin = .data$.lower, ymax = .data$.upper,
            fill = if (by) .data$.by), colour = NA, alpha = 0.2,
            na.rm = TRUE) +
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
            colour = group.by, fill = group.by) +
        .scale_categorical(c("colour", "fill")) + .theme_tempoPlot()
}
