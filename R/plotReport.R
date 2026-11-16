#' Plot a report of aligned panels
#'
#' \code{plotReport()} stacks tempoPlot panels over a shared time axis for
#' one subject, a cohort or an intervention.
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
#' @param y \code{Character vector}. Variables drawn as trends, one panel
#'   each: columns of \code{colData(x)}, features of \code{assay.type} or
#'   \code{"<dimred>:<component>"}. (Default: \code{NULL})
#'
#' @param level \code{Character scalar}. Draws one subject (\code{"subject"}):
#'   its trends and composition; a cohort (\code{"cohort"}): the sampling
#'   timeline, boxes of the trends at each time and the mean composition;
#'   or an intervention (\code{"intervention"}): the change of the trends
#'   from each subject's first sample with the mean of each
#'   \code{group.by} group, and the mean composition of each group.
#'   (Default: \code{"subject"})
#'
#' @param subject \code{Character scalar}. The subject drawn with
#'   \code{level = "subject"}. (Default: \code{NULL})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} with the arm of each subject, needed for
#'   \code{level = "intervention"}. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of
#'   abundances for the composition, and of features in \code{y}.
#'   (Default: \code{NULL}, the first assay)
#'
#' @param rank \code{Character scalar}. Specifies the column of
#'   \code{rowData(x)} whose levels are summed in the composition.
#'   (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features shown in the
#'   composition; the rest are summed into "Other". (Default: \code{5})
#'
#' @param ... additional arguments passed to
#'   \code{\link[patchwork:wrap_plots]{patchwork::wrap_plots}}.
#'
#' @details
#' Requires patchwork. The panels are \code{\link{plotSwimmer}},
#' \code{\link{plotTrend}} and \code{\link{plotComposition}}; build other
#' combinations by adding their plots with patchwork.
#'
#' @return
#' A \code{patchwork} object.
#'
#' @name plotReport
#'
#' @examples
#' # One bioreactor: short-chain fatty acids and species
#' data("minimalgut", package = "miaTime")
#' plotReport(minimalgut, time.col = "Time_hr_num",
#'     group = "StudyIdentifier", y = c("Acetate", "Butyrate"),
#'     subject = "Bioreactor A", top = 8)
#'
#' # All bioreactors
#' plotReport(minimalgut, time.col = "Time_hr_num",
#'     group = "StudyIdentifier", y = "Acetate", level = "cohort", top = 8)
NULL

#' @rdname plotReport
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotReport", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, y = NULL,
            level = c("subject", "cohort", "intervention"), subject = NULL,
            group.by = NULL, assay.type = NULL, rank = NULL, top = 5, ...) {
        level <- match.arg(level)
        if (!requireNamespace("patchwork", quietly = TRUE)) {
            stop("plotReport() requires the patchwork package.",
                call. = FALSE)
        }
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_column(x, group)
        if (!is.null(y) && !is.character(y)) {
            stop("'y' must be NULL or a character vector.", call. = FALSE)
        }
        if (level == "subject") {
            if (!.is_one_name(subject) ||
                    !subject %in% colData(x)[[group]]) {
                stop("'subject' must be one value of the 'group' column.",
                    call. = FALSE)
            }
            x <- x[, colData(x)[[group]] %in% subject]
        }
        if (level == "intervention" && is.null(group.by)) {
            stop("level = \"intervention\" needs 'group.by'.",
                call. = FALSE)
        }
        trend <- function(v) {
            switch(level,
                subject = plotTrend(x, v, time.col, group,
                    assay.type = assay.type),
                cohort = plotTrend(x, v, time.col, group,
                    assay.type = assay.type, smooth = "mean"),
                intervention = plotTrend(x, v, time.col, group,
                    assay.type = assay.type, baseline = "difference",
                    smooth = "mean", group.by = group.by))
        }
        panels <- lapply(y, trend)
        if (level == "cohort") {
            panels <- c(list(plotSwimmer(x, time.col, group)), panels)
        }
        if (level == "cohort") {
            x$.all <- "all subjects"
        }
        panels <- c(panels, list(plotComposition(x, time.col,
            group = group,
            assay.type = assay.type, rank = rank, top = top,
            layout = "area",
            level = if (level == "subject") "subject" else "group",
            group.by = switch(level, subject = NULL, cohort = ".all",
                intervention = group.by))))
        patchwork::wrap_plots(panels, ncol = 1L, ...)
    }
)
