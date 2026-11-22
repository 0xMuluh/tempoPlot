#' Plot an effect over time
#'
#' \code{plotTimeEffect()} draws a model estimate at each time with its
#' interval, against a reference line at zero: a difference between groups
#' at each time point, a time by group interaction, or a GAM difference
#' curve.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object with the results in \code{metadata(x)}, or a \code{data.frame}
#'   of results with the columns \code{time.col} and \code{estimate}, and
#'   optionally \code{lower} and \code{upper} for an interval.
#'
#' @param name \code{Character scalar}. Specifies the entry of
#'   \code{metadata(x)} holding the results.
#'
#' @param time.col \code{Character scalar}. Specifies the column of the
#'   results that contains the times: numeric, \code{Date} or a factor of
#'   visits.
#'
#' @param group.by \code{Character scalar}. Specifies the column of the
#'   results that splits the estimates into curves, such as the contrast or
#'   the feature. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' The y axis shows the effect, not the value of a variable; to draw fitted
#' values over the data, use the \code{fit} argument of
#' \code{\link{plotTrend}}.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotTimeEffect
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#'
#' # Difference in diversity between men and women at each visit, from a
#' # linear model per visit; each subject has one sample per visit
#' cd <- as.data.frame(SummarizedExperiment::colData(tse))
#' effect <- do.call(rbind, lapply(split(cd, cd$visit), function(d) {
#'     model <- lm(diversity ~ sex, data = d)
#'     ci <- confint(model)["sexmale", ]
#'     data.frame(visit = d$visit[1L], estimate = coef(model)[["sexmale"]],
#'         lower = ci[[1L]], upper = ci[[2L]])
#' }))
#' plotTimeEffect(effect, time.col = "visit")
#'
#' # The same results stored in metadata(tse)
#' S4Vectors::metadata(tse)$sex_effect <- effect
#' plotTimeEffect(tse, name = "sex_effect", time.col = "visit")
NULL

#' @rdname plotTimeEffect
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotTimeEffect", signature = c(x = "SummarizedExperiment"),
    function(x, name, time.col, group.by = NULL, ...) {
        if (!.is_one_name(name) || !is.data.frame(metadata(x)[[name]])) {
            stop("'name' must name a table in metadata(x).", call. = FALSE)
        }
        plotTimeEffect(metadata(x)[[name]], time.col, group.by, ...)
    }
)

#' @rdname plotTimeEffect
#' @importFrom methods setMethod
#' @export
setMethod("plotTimeEffect", signature = c(x = "data.frame"),
    function(x, time.col, group.by = NULL, ...) {
        if (!.is_one_name(time.col) || !time.col %in% names(x)) {
            stop("'time.col' must be one column name of 'x'.", call. = FALSE)
        }
        if (!is.null(group.by) &&
                !(.is_one_name(group.by) && group.by %in% names(x))) {
            stop("'group.by' must be NULL or one column name of 'x'.",
                call. = FALSE)
        }
        time <- x[[time.col]]
        if (!(is.numeric(time) || inherits(time, "Date") || is.factor(time))) {
            stop("'time.col' must be numeric, Date or factor.", call. = FALSE)
        }
        if (!is.numeric(x$estimate) ||
                sum(c("lower", "upper") %in% names(x)) == 1L) {
            stop("'x' must have a numeric column estimate, and lower and ",
                "upper together.", call. = FALSE)
        }
        .plot_time_effect(as.data.frame(x), time.col, group.by, ...)
    }
)

#' @importFrom ggplot2 ggplot geom_hline labs
.plot_time_effect <- function(df, time.col, group.by, ...) {
    ggplot(df) +
        geom_hline(yintercept = 0, colour = "grey40") +
        .fit_layers(df, time.col, group.by, ...) +
        labs(x = time.col, y = "estimate", colour = group.by)
}
