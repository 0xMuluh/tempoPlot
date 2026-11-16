#' Plot a variogram against time lag
#'
#' \code{plotLagProfile()} draws how far samples of the same subject are
#' apart at each time lag: the variogram \eqn{D(L)} with its confidence
#' interval, as stored by \code{miaTime::addVariogram()}.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param name \code{Character scalar}. Specifies the entry of
#'   \code{metadata(x)} holding the variogram. (Default: \code{"variogram"})
#'
#' @param colour.by \code{Character scalar}. Specifies the column of the
#'   variogram that colours the estimates, such as the \code{split.by}
#'   variable of \code{addVariogram()}. (Default: \code{NULL})
#'
#' @param type \code{Character scalar}. Draws the dissimilarity at each lag
#'   (\code{"dissimilarity"}) or, for a variogram computed with
#'   \code{split.by}, the ratio between its levels (\code{"ratio"}), on a
#'   log scale with a line at 1. (Default: \code{"dissimilarity"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_linerange]{geom_pointrange}}.
#'
#' @details
#' Lags the design does not support (missing estimates) are dropped. With
#' \code{split.by}, the estimates of the levels are slightly shifted apart
#' at each lag.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotLagProfile
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- mia::transformAssay(hitchip1006, method = "hellinger")
#'
#' # Squared Hellinger distance between samples of a subject L time units
#' # apart
#' tse <- miaTime::addVariogram(tse, time.col = "time", group = "subject",
#'     assay.type = "hellinger", method = "euclidean", squared = TRUE,
#'     time.lag = c(1.5, 2, 2.5, 3, 4, 6, 8, 9))
#' plotLagProfile(tse)
#'
#' # The same for women and men
#' tse <- miaTime::addVariogram(tse, time.col = "time", group = "subject",
#'     assay.type = "hellinger", method = "euclidean", squared = TRUE,
#'     time.lag = c(2, 6, 8), split.by = "sex", name = "variogram_sex")
#' plotLagProfile(tse, name = "variogram_sex", colour.by = "sex")
#'
#' # Their ratio
#' plotLagProfile(tse, name = "variogram_sex", type = "ratio")
NULL

#' @rdname plotLagProfile
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotLagProfile", signature = c(x = "SummarizedExperiment"),
    function(x, name = "variogram", colour.by = NULL,
            type = c("dissimilarity", "ratio"), ...) {
        type <- match.arg(type)
        if (!.is_one_name(name) || is.null(metadata(x)[[name]])) {
            stop("'name' must be an entry of metadata(x); compute it with ",
                "miaTime::addVariogram().", call. = FALSE)
        }
        df <- as.data.frame(metadata(x)[[name]])
        if (!all(c("lag", "estimate", "lower", "upper") %in% names(df))) {
            stop("metadata(x)$", name, " must have the columns lag, ",
                "estimate, lower and upper.", call. = FALSE)
        }
        if (!is.null(colour.by) && !(.is_one_name(colour.by) &&
                colour.by %in% names(df))) {
            stop("'colour.by' must be NULL or one column of metadata(x)$",
                name, ".", call. = FALSE)
        }
        if ("type" %in% names(df)) {
            df <- df[df$type == type, , drop = FALSE]
        } else if (type == "ratio") {
            df <- df[0L, , drop = FALSE]
        }
        df <- df[!is.na(df$estimate), , drop = FALSE]
        if (nrow(df) == 0L) {
            stop("metadata(x)$", name, " has no ", type, " estimates",
                if (type == "ratio") "; compute it with 'split.by'", ".",
                call. = FALSE)
        }

        # The split.by column sits between lag and type; its ratio rows
        # name the levels, such as "male / female".
        ylab <- "D(L)"
        if (type == "ratio") {
            split <- setdiff(names(df)[seq_len(match("type", names(df)))],
                c("lag", "type"))
            ylab <- paste("Ratio of D(L)", if (length(split) == 1L)
                paste0("(", unique(df[[split]])[1L], ")"))
        }
        .plot_lag_profile(df, colour.by, type, ylab, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_pointrange labs
#'   position_dodge geom_hline scale_y_log10
.plot_lag_profile <- function(df, colour.by, type = "dissimilarity",
        ylab = "D(L)", ...) {
    # Groups at the same lag are shifted apart so their intervals show.
    dodge <- position_dodge(width = 0.02 * diff(range(df$lag)))
    p <- ggplot(df, aes(x = .data$lag, y = .data$estimate,
            colour = if (!is.null(colour.by)) .data[[colour.by]],
            group = if (!is.null(colour.by)) .data[[colour.by]] else 1L)) +
        geom_line(position = dodge) +
        geom_pointrange(aes(ymin = .data$lower, ymax = .data$upper),
            position = dodge, ...) +
        labs(x = "lag", y = ylab, colour = colour.by)
    if (type == "ratio") {
        # A ratio of 1 means the same dissimilarity in both levels.
        p <- p + geom_hline(yintercept = 1, linetype = "dashed",
                colour = "grey40") +
            scale_y_log10()
    }
    p
}
