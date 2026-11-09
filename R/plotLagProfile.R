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
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_linerange]{geom_pointrange}}.
#'
#' @details
#' Lags the design does not support (missing estimates) are dropped. With
#' \code{split.by}, only the estimates of each level are drawn, not the
#' ratios between levels, slightly shifted apart at each lag.
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
NULL

#' @rdname plotLagProfile
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotLagProfile", signature = c(x = "SummarizedExperiment"),
    function(x, name = "variogram", colour.by = NULL, ...) {
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
            df <- df[df$type == "dissimilarity", , drop = FALSE]
        }
        df <- df[!is.na(df$estimate), , drop = FALSE]

        .plot_lag_profile(df, colour.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_pointrange labs
#'   position_dodge
.plot_lag_profile <- function(df, colour.by, ...) {
    # Groups at the same lag are shifted apart so their intervals show.
    dodge <- position_dodge(width = 0.02 * diff(range(df$lag)))
    ggplot(df, aes(x = .data$lag, y = .data$estimate,
            colour = if (!is.null(colour.by)) .data[[colour.by]],
            group = if (!is.null(colour.by)) .data[[colour.by]] else 1L)) +
        geom_line(position = dodge) +
        geom_pointrange(aes(ymin = .data$lower, ymax = .data$upper),
            position = dodge, ...) +
        labs(x = "lag", y = "D(L)", colour = colour.by)
}
