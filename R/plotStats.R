#' Plot ranked statistics of features or subjects
#'
#' \code{plotStats()} draws one statistic per feature or per subject, ranked
#' by value, with optional intervals. It plots results stored by other
#' functions, such as bimodality and stability scores or Cox coefficients
#' from miaTime, and per-subject estimates such as slopes.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param stat \code{Character scalar}. Specifies the column holding the
#'   statistic: a column of \code{rowData(x)} for features, or of
#'   \code{colData(x)} for subjects when \code{group} is given.
#'
#' @param lower,upper \code{Character scalar}. Specify the columns holding
#'   the interval around \code{stat}. Give both or neither.
#'   (Default: \code{NULL})
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If given, one row is drawn
#'   per subject, and \code{stat}, \code{lower}, \code{upper} and
#'   \code{colour.by} must be constant within each subject.
#'   (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Keeps the features or subjects with the
#'   largest absolute statistics. If \code{NULL}, keeps all.
#'   (Default: \code{NULL})
#'
#' @param colour.by \code{Character scalar}. Specifies the column, from the
#'   same table as \code{stat}, that colours the points.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Rows are ordered by the statistic, the largest at the top. Missing
#' statistics are dropped.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotStats
#'
#' @examples
#' # Bimodality of each genus over all samples, computed by miaTime
#' data("hitchip1006", package = "miaTime")
#' tse <- miaTime::addBimodality(hitchip1006, assay.type = "counts")
#' plotStats(tse, stat = "bimodality", top = 20, colour.by = "Phylum")
#'
#' # Cox coefficients for time to relapse, computed by miaTime
#' data("crohn_survival", package = "miaTime")
#' set.seed(1)
#' tse <- miaTime::addSurvival(crohn_survival, assay.type = "counts",
#'     time.col = "event_time", event.col = "event")
#' plotStats(tse, stat = "survival_coef", top = 10)
#'
#' # Mean diversity of each subject with its range over visits
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
#' tse$mean <- ave(tse$diversity, tse$subject, FUN = mean)
#' tse$min <- ave(tse$diversity, tse$subject, FUN = min)
#' tse$max <- ave(tse$diversity, tse$subject, FUN = max)
#' plotStats(tse, stat = "mean", lower = "min", upper = "max",
#'     group = "subject", colour.by = "sex")
NULL

#' @rdname plotStats
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom methods setMethod
#' @export
setMethod("plotStats", signature = c(x = "SummarizedExperiment"),
    function(x, stat, lower = NULL, upper = NULL, group = NULL, top = NULL,
            colour.by = NULL, ...) {
        where <- if (is.null(group)) "rowData" else "colData"
        if (is.null(stat)) {
            stop("'stat' must be one column name of ", where, "(x).",
                call. = FALSE)
        }
        .check_column(x, group)
        .check_column(x, stat, where = where)
        .check_column(x, lower, where = where)
        .check_column(x, upper, where = where)
        .check_column(x, colour.by, where = where)
        if (xor(is.null(lower), is.null(upper))) {
            stop("'lower' and 'upper' must be given together.",
                call. = FALSE)
        }
        if (!is.null(top) && !(is.numeric(top) && length(top) == 1L &&
                top >= 1)) {
            stop("'top' must be NULL or a single number of at least 1.",
                call. = FALSE)
        }

        cols <- unique(c(group, stat, lower, upper, colour.by))
        tab <- if (is.null(group)) rowData(x) else colData(x)
        df <- as.data.frame(tab[, cols, drop = FALSE], optional = TRUE)
        df$.id <- if (is.null(group)) rownames(x) else df[[group]]
        df <- unique(df)
        if (anyDuplicated(df$.id)) {
            stop("'stat', 'lower', 'upper' and 'colour.by' must be constant ",
                "within each subject of 'group'.", call. = FALSE)
        }
        if (!all(vapply(df[c(stat, lower, upper)], is.numeric, logical(1)))) {
            stop("'stat', 'lower' and 'upper' must be numeric columns.",
                call. = FALSE)
        }

        df <- df[!is.na(df[[stat]]), , drop = FALSE]
        df <- df[order(-abs(df[[stat]])), , drop = FALSE]
        df <- df[seq_len(min(top %||% nrow(df), nrow(df))), , drop = FALSE]
        df$.id <- factor(df$.id, levels = df$.id[order(df[[stat]])])

        .plot_stats(df, stat, lower, upper, group, colour.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_point geom_linerange labs
.plot_stats <- function(df, stat, lower, upper, group, colour.by, ...) {
    p <- ggplot(df, aes(x = .data[[stat]], y = .data$.id,
        colour = if (!is.null(colour.by)) .data[[colour.by]]))
    if (!is.null(lower)) {
        p <- p + geom_linerange(aes(xmin = .data[[lower]],
            xmax = .data[[upper]]))
    }
    p + geom_point(...) + labs(x = stat, y = group, colour = colour.by)
}
