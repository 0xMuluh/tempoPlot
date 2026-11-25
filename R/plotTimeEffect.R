#' Plot an effect over time
#'
#' \code{plotTimeEffect()} draws a model estimate at each time with its
#' interval, against a line at its value without an effect (zero by
#' default): a difference between groups at each time point, a time by
#' group interaction, or a GAM difference curve.
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
#'   results that contains the times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits.
#'
#' @param group.by \code{Character scalar}. Specifies the column of the
#'   results that splits the estimates into curves, such as the contrast or
#'   the feature. (Default: \code{NULL})
#'
#' @param null \code{data.frame}. Effects estimated under the null
#'   hypothesis, such as after permuting the group labels, with the
#'   columns \code{permutation}, \code{time.col} and \code{estimate}, and
#'   optionally \code{group.by}. Each permutation is drawn as a thin
#'   translucent line behind the observed effect. (Default: \code{NULL})
#'
#' @param null.value \code{Numeric scalar}. The value of the effect
#'   without an effect, drawn as a line, such as 1 for a ratio.
#'   (Default: \code{0})
#'
#' @param p.col \code{Character scalar}. Specifies the column of p-values
#'   to draw on a log scale instead of the effect, with a dashed line at
#'   \code{threshold}. (Default: \code{NULL})
#'
#' @param threshold \code{Numeric scalar}. Significance level for
#'   \code{p.col}. (Default: \code{0.05})
#'
#' @param size.by \code{Character scalar}. Specifies the column that sets
#'   the size of a point at each time, such as the number of subjects.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' The y axis shows the effect, not the value of a variable; to draw fitted
#' values over the data, use the \code{fit} argument of
#' \code{\link{plotTrend}}.
#'
#' The \code{null} curves and the p-value curve with points sized by the
#' number of subjects follow the plots of SplinectomeR (Shields-Cutler et
#' al. 2018, doi:10.3389/fmicb.2018.00785, Fig. 2).
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
#' sex_effect <- function(cd) {
#'     do.call(rbind, lapply(split(cd, cd$visit), function(d) {
#'         model <- lm(diversity ~ sex, data = d)
#'         ci <- confint(model)["sexmale", ]
#'         data.frame(visit = d$visit[1L],
#'             estimate = coef(model)[["sexmale"]],
#'             lower = ci[[1L]], upper = ci[[2L]],
#'             p = summary(model)$coefficients["sexmale", 4L], n = nrow(d))
#'     }))
#' }
#' cd <- as.data.frame(SummarizedExperiment::colData(tse))
#' effect <- sex_effect(cd)
#' plotTimeEffect(effect, time.col = "visit")
#'
#' # Against the effects after permuting sex within each visit 50 times
#' set.seed(1)
#' null <- do.call(rbind, lapply(seq_len(50L), function(i) {
#'     cd$sex <- ave(as.character(cd$sex), cd$visit, FUN = sample)
#'     data.frame(permutation = i, sex_effect(cd))
#' }))
#' plotTimeEffect(effect, time.col = "visit", null = null)
#'
#' # p-values, with points sized by the number of subjects
#' plotTimeEffect(effect, time.col = "visit", p.col = "p", size.by = "n")
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
    function(x, time.col, group.by = NULL, null = NULL, null.value = 0,
            p.col = NULL, threshold = 0.05, size.by = NULL, ...) {
        .check_number(null.value)
        .check_table_columns(x, time.col = time.col, null = FALSE)
        .check_table_columns(x, group.by = group.by, p.col = p.col,
            size.by = size.by)
        type <- .time_type(x[[time.col]])
        if (is.null(p.col)) {
            if (!is.numeric(x$estimate) ||
                    sum(c("lower", "upper") %in% names(x)) == 1L) {
                stop("'x' must have a numeric column estimate, and lower ",
                    "and upper together.", call. = FALSE)
            }
        } else {
            .check_number(threshold, 0, 1, open = TRUE)
            if (!is.null(null)) {
                stop("'null' must be NULL when 'p.col' is set.",
                    call. = FALSE)
            }
            # The p-values take the place of the effect.
            x$estimate <- x[[p.col]]
            x$lower <- x$upper <- NULL
        }
        if (!is.null(null) && !(is.data.frame(null) &&
                all(c("permutation", time.col) %in% names(null)) &&
                is.numeric(null$estimate) &&
                identical(.time_type(null[[time.col]]), type))) {
            stop("'null' must be NULL or a data.frame with columns ",
                "permutation, '", time.col, "' of the type of 'x', and a ",
                "numeric estimate.", call. = FALSE)
        }
        .plot_time_effect(as.data.frame(x), time.col, group.by, null,
            null.value, p.col, threshold, size.by, ...)
    }
)

#' @importFrom ggplot2 ggplot geom_hline geom_line geom_point aes .data
#'   labs scale_y_log10
.plot_time_effect <- function(df, time.col, group.by, null, null.value,
        p.col, threshold, size.by, ...) {
    p <- ggplot(df) + if (is.null(p.col)) {
        geom_hline(yintercept = null.value, colour = "grey40")
    } else {
        list(geom_hline(yintercept = threshold, colour = "grey40",
            linetype = "dashed"), scale_y_log10())
    }
    if (!is.null(null)) {
        p <- p + .null_layer(as.data.frame(null), time.col, group.by)
    }
    p <- p + .fit_layers(df, time.col, group.by, ...)
    if (!is.null(size.by)) {
        p <- p + geom_point(aes(x = .data[[time.col]], y = .data$estimate,
            size = .data[[size.by]],
            colour = if (!is.null(group.by)) .data[[group.by]]))
    }
    p + labs(x = time.col, y = if (is.null(p.col)) "estimate" else p.col,
            colour = group.by, size = size.by) +
        .scale_by(df, group.by, c("colour", "fill")) +
        .theme_tempoPlot()
}

# One thin translucent line per permutation of null (and group.by, if
# null has it, in the colour of its group; else grey).
.null_layer <- function(null, time.col, group.by) {
    by <- if (!is.null(group.by) && group.by %in% names(null)) group.by
    null$.curve <- paste(null$permutation, if (!is.null(by)) null[[by]])
    do.call(geom_line, c(list(aes(x = .data[[time.col]], y = .data$estimate,
            group = .data$.curve, colour = if (!is.null(by)) .data[[by]]),
        data = null, inherit.aes = FALSE, alpha = 0.2, linewidth = 0.3,
        show.legend = FALSE), if (is.null(by)) list(colour = .na_colour)))
}
