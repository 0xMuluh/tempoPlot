#' Plot two subject-level statistics against each other
#'
#' \code{plotSubjectScatter()} draws one point per subject from two
#' statistics that are constant within each subject, such as mean diversity
#' and its variability, resistance and resilience, or subject loadings of a
#' decomposition stored in \code{colData(x)}.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param x.stat,y.stat \code{Character scalar}. Specify the numeric columns
#'   from \code{colData(x)} drawn on the x and y axes.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param colour.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)}, constant within subjects, that colours the points.
#'   (Default: \code{NULL})
#'
#' @param quadrants \code{Logical scalar} or \code{Numeric vector} of
#'   length 2. Draws reference lines that split the plot into quadrants: at
#'   the medians of the two statistics (\code{TRUE}) or at the given x and y
#'   values. (Default: \code{FALSE})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Subjects with a missing statistic are dropped.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotSubjectScatter
#'
#' @examples
#' # Mean diversity against its standard deviation over each subject's
#' # samples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
#' tse$mean <- ave(tse$diversity, tse$subject, FUN = mean)
#' tse$sd <- ave(tse$diversity, tse$subject, FUN = sd)
#' plotSubjectScatter(tse, x.stat = "mean", y.stat = "sd",
#'     group = "subject", colour.by = "sex", quadrants = TRUE)
NULL

#' @rdname plotSubjectScatter
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median
#' @importFrom methods setMethod
#' @export
setMethod("plotSubjectScatter", signature = c(x = "SummarizedExperiment"),
    function(x, x.stat, y.stat, group, colour.by = NULL, quadrants = FALSE,
            ...) {
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        for (arg in c("x.stat", "y.stat")) {
            col <- get(arg)
            if (is.null(col)) {
                stop("'", arg, "' must be one column name of colData(x).",
                    call. = FALSE)
            }
            .check_column(x, col, arg = arg)
            if (!is.numeric(colData(x)[[col]])) {
                stop("'", arg, "' must be a numeric column.", call. = FALSE)
            }
        }
        .check_column(x, group)
        .check_column(x, colour.by)
        if (!(isTRUE(quadrants) || isFALSE(quadrants) ||
                (is.numeric(quadrants) && length(quadrants) == 2L))) {
            stop("'quadrants' must be TRUE, FALSE or two numbers.",
                call. = FALSE)
        }

        cols <- unique(c(group, x.stat, y.stat, colour.by))
        df <- unique(as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE))
        if (anyDuplicated(df[[group]])) {
            stop("'x.stat', 'y.stat' and 'colour.by' must be constant ",
                "within each subject of 'group'.", call. = FALSE)
        }
        df <- df[!is.na(df[[x.stat]]) & !is.na(df[[y.stat]]), ,
            drop = FALSE]
        if (nrow(df) == 0L) {
            stop("No subject has both statistics.", call. = FALSE)
        }
        if (isTRUE(quadrants)) {
            quadrants <- c(median(df[[x.stat]]), median(df[[y.stat]]))
        }
        .plot_subject_scatter(df, x.stat, y.stat, colour.by, quadrants, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_point geom_hline geom_vline
#'   labs
.plot_subject_scatter <- function(df, x.stat, y.stat, colour.by, quadrants,
        ...) {
    p <- ggplot(df, aes(x = .data[[x.stat]], y = .data[[y.stat]]))
    if (is.numeric(quadrants)) {
        p <- p + geom_vline(xintercept = quadrants[[1L]], linetype = "dashed",
                colour = "grey40") +
            geom_hline(yintercept = quadrants[[2L]], linetype = "dashed",
                colour = "grey40")
    }
    p + geom_point(aes(colour = if (!is.null(colour.by))
            .data[[colour.by]]), ...) +
        labs(x = x.stat, y = y.stat, colour = colour.by)
}
