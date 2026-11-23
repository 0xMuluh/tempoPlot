#' Plot how long features stay present or subjects stay in a state
#'
#' \code{plotRuns()} finds the runs of each subject, consecutive samples in
#' which a feature stays present or absent, or a state holds, and draws
#' the distribution of their durations as Kaplan-Meier curves, one panel
#' per present, absent or state.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param state.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} with the state of each sample, such as a community
#'   type. If \code{NULL}, runs are of the presence and absence of
#'   features. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay for the
#'   presence of features. (Default: \code{"counts"})
#'
#' @param features \code{Character vector}. Features whose runs are pooled.
#'   If \code{NULL}, pools all features. (Default: \code{NULL})
#'
#' @param detection \code{Numeric scalar}. Specifies the value above which
#'   a feature is present. (Default: \code{0})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the runs into curves, taken from the
#'   first sample of each run. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Ends a run where consecutive
#'   samples of a subject are further apart than this, in \code{time.col}
#'   units. For visits, \code{NULL} ends runs at visits without a sample;
#'   otherwise \code{NULL} never ends them. (Default: \code{NULL})
#'
#' @param layout \code{Character scalar}. Draws the Kaplan-Meier curves
#'   with confidence bands and censoring marks (\code{"survival"}) or the
#'   same curves on log-log axes, the complementary cumulative
#'   distribution (\code{"ccdf"}). (Default: \code{"survival"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_step]{geom_step}}.
#'
#' @details
#' A run lasts from its first sample to the first sample of the next run,
#' the first time the change is seen, as in \code{\link{plotSurvival}}.
#' A run that reaches the last sample of a subject, or a gap longer than
#' \code{max.gap}, has no observed end: it is censored at its last sample
#' and marked with a cross. A run that starts at the first sample of a
#' subject, or after such a gap, may have begun earlier; its duration
#' counts from that sample. Samples without a state, or with a missing
#' assay value, are skipped. Durations are in \code{time.col} units: days
#' for \code{Date}, hours for \code{POSIXct} times and visit steps for
#' visits.
#'
#' Runs of presence give the residence times of features and runs of
#' absence their return times, pooled over \code{features}. With
#' \code{layout = "ccdf"}, the curves are drawn on log-log axes, as the
#' residence and return time distributions of Ji et al. (2020); runs of
#' length zero and the drop to zero after the longest run are not drawn.
#' Runs of a state give the durations of its episodes, as the dysbiotic
#' episodes and the intervals between them of Lloyd-Price et al. (2019).
#'
#' @references
#' Ji BW, et al. (2020). Macroecological dynamics of gut microbiota.
#' \emph{Nature Microbiology} 5(5):768-775.
#' doi:10.1038/s41564-020-0685-1
#'
#' Lloyd-Price J, et al. (2019). Multi-omics of the gut microbial ecosystem
#' in inflammatory bowel diseases. \emph{Nature} 569(7758):655-662.
#' doi:10.1038/s41586-019-1237-9
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotRuns
#'
#' @examples
#' # Residence and return times of the species in three bioreactors
#' data("minimalgut", package = "miaTime")
#' plotRuns(minimalgut, time.col = "Time.hr", group = "StudyIdentifier",
#'     layout = "ccdf")
#'
#' # Days in each enterotype, ending runs at gaps of more than 3 days.
#' # Samples without a cell count have the enterotype "NA".
#' data("temporalMicrobiome20", package = "miaTime")
#' tse <- temporalMicrobiome20[, !is.na(temporalMicrobiome20$day)]
#' tse$enterotype[tse$enterotype == "NA"] <- NA
#' plotRuns(tse, time.col = "day", group = "person_id",
#'     state.col = "enterotype", max.gap = 3)
NULL

#' @rdname plotRuns
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotRuns", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, state.col = NULL, assay.type = "counts",
            features = NULL, detection = 0, group.by = NULL,
            max.gap = NULL, layout = c("survival", "ccdf"), ...) {
        layout <- match.arg(layout)
        .require_column(x, group)
        .check_longitudinal(x, time.col, group)
        .check_column(x, group.by)
        .check_column(x, state.col)

        if (is.null(state.col)) {
            mat <- .get_assay(x, assay.type)
            if (!is.null(features)) {
                mat <- mat[.select_features(mat, features, NULL), ,
                    drop = FALSE]
            }
            present <- .present(mat, detection)
            values <- lapply(seq_len(nrow(present)), function(j)
                c("absent", "present")[present[j, ] + 1L])
            states <- c("present", "absent")
        } else {
            state <- colData(x)[[state.col]]
            values <- list(as.character(state))
            states <- .levels_of(state)
        }

        df <- as.data.frame(colData(x)[, c(time.col, group), drop = FALSE],
            optional = TRUE)
        df$.i <- seq_len(ncol(x))
        df <- .order_time(df, time.col, group)
        df$.t <- .time_numeric(df[[time.col]])
        visit <- is.factor(df[[time.col]])
        runs <- do.call(rbind, lapply(values, function(v)
            .find_runs(v[df$.i], df, group, max.gap, visit)))
        if (is.null(runs) || nrow(runs) == 0L) {
            stop("No sample has a state or an assay value.", call. = FALSE)
        }
        runs$.state <- factor(runs$.state, levels = states)
        runs$.by <- if (is.null(group.by)) "" else
            as.character(colData(x)[[group.by]][runs$.i])

        km <- do.call(rbind, lapply(split(runs, runs$.state, drop = TRUE),
            function(d) cbind(.kaplan_meier(d), .state = d$.state[1L])))
        .plot_runs(km, time.col, group.by, layout, ...)
    }
)

# Runs of equal values v in the samples of df, sorted by .order_time(),
# with numeric times .t. NA values are skipped; a run also ends at a new
# subject or a gap above max.gap. One row per run: its value (.state),
# duration (.time), whether its end was seen (.status) and the column of
# x of its first sample (.i).
.find_runs <- function(v, df, group, max.gap, visit) {
    keep <- !is.na(v)
    v <- v[keep]
    df <- .break_gaps(df[keep, , drop = FALSE], ".t", group, max.gap,
        visit = visit)
    n <- length(v)
    # link: sample k and k + 1 are consecutive samples of one series.
    link <- seq_len(n) %in% .next_in_subject(df[[group]],
        segment = df$.segment)
    last <- which(!(link & c(v[-1L] == v[-n], FALSE)))
    first <- c(1L, last + 1L)[seq_along(last)]
    seen <- link[last]
    data.frame(.state = v[first], .time = df$.t[last + seen] - df$.t[first],
        .status = as.integer(seen), .i = df$.i[first])
}

#' @importFrom ggplot2 ggplot aes .data geom_step labs facet_wrap vars
#'   scale_x_log10 scale_y_log10
.plot_runs <- function(km, time.col, group.by, layout, ...) {
    x.title <- paste(time.col, "from the start of the run")
    y.title <- "share of runs lasting longer"
    if (layout == "survival") {
        p <- .plot_survival(km, group.by, x.title, y.title, ...)
    } else {
        km <- km[km$.time > 0 & km$.surv > 0, , drop = FALSE]
        p <- ggplot(km, aes(x = .data$.time, y = .data$.surv,
                colour = if (!is.null(group.by)) .data$.by)) +
            geom_step(...) +
            scale_x_log10() + scale_y_log10() +
            labs(x = x.title, y = y.title, colour = group.by) +
            .scale_categorical("colour") + .theme_tempoPlot()
    }
    p + facet_wrap(vars(.data$.state))
}
