#' Plot the distribution of abundance changes between samples
#'
#' \code{plotChanges()} pools the log ratios
#' \eqn{\log(x_{t+1} / x_t)} of each feature between consecutive samples
#' of a subject and draws their distribution on a log density axis.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits. Time differences, as in
#'   \code{max.gap}, are in days for \code{Date} and hours for
#'   \code{POSIXct} times.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay.
#'   (Default: \code{"counts"})
#'
#' @param detection \code{Numeric scalar}. A change is counted when the
#'   feature is above this in both samples. (Default: \code{0})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the changes into one distribution per
#'   group, by its value in the earlier sample. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Skips the intervals between
#'   consecutive samples of a subject that are further apart than this, in
#'   \code{time.col} units. For visits, \code{NULL} skips intervals across
#'   visits without a sample; otherwise \code{NULL} keeps all of them.
#'   (Default: \code{NULL})
#'
#' @param bins \code{Integer scalar}. Approximate number of bins of the
#'   log ratios. (Default: \code{30})
#'
#' @param fit \code{Character scalar}. Draws no fit (\code{"none"}) or the
#'   Laplace distribution fitted by maximum likelihood (\code{"laplace"}),
#'   as a dashed line. (Default: \code{"none"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Ji et al. (2020, doi:10.1038/s41564-020-0685-1, Fig. 1) found that daily
#' changes of gut bacteria follow a Laplace distribution, whose log density
#' falls linearly on both sides; a normal distribution would fall as a
#' parabola. The Laplace fit has the median as location and the mean
#' absolute deviation from it as scale \eqn{b}, given in the subtitle.
#' Each point is the density of a bin with at least one change; empty bins
#' are not drawn.
#'
#' Changes over longer intervals tend to be larger, so a distribution is
#' comparable to another only at the same interval. The subtitle gives the
#' range of intervals; use \code{max.gap} to keep only the shortest
#' intervals of regularly sampled data, such as \code{max.gap = 1} for
#' daily samples. The ratios of counts include the change in sequencing
#' depth; use relative or absolute abundances instead.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotChanges
#'
#' @examples
#' # Daily changes of the absolute abundances of genera in 20 women, with
#' # the Laplace distribution of the same median and mean absolute
#' # deviation
#' data("temporalMicrobiome20", package = "miaTime")
#' tse <- temporalMicrobiome20
#' # Three samples have no day
#' tse <- tse[, !is.na(tse$day)]
#' plotChanges(tse, time.col = "day", group = "person_id", max.gap = 1,
#'     fit = "laplace")
#'
#' # Changes of relative abundances over at most 4 hours in each of three
#' # bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "relabundance")
#' plotChanges(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     assay.type = "relabundance", max.gap = 4,
#'     group.by = "StudyIdentifier")
NULL

#' @rdname plotChanges
#' @importFrom methods setMethod
#' @export
setMethod("plotChanges", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            detection = 0, group.by = NULL, max.gap = NULL, bins = 30,
            fit = c("none", "laplace"), ...) {
        fit <- match.arg(fit)
        .check_number(bins, 1)
        df <- .change_table(x, time.col, group, assay.type, detection,
            group.by, max.gap)
        breaks <- pretty(range(df$.change), n = bins)
        res <- lapply(split(df$.change, df$.by), function(v) {
            list(bins = .change_density(v, breaks),
                curve = if (fit == "laplace") .fit_laplace(v))
        })
        .plot_changes(df, .bind_groups(res, "bins"),
            .bind_groups(res, "curve"), time.col, group.by, ...)
    }
)

# One row per feature and pair of consecutive samples of a subject where
# the feature is above detection in both: the log ratio .change, the
# interval .interval and the group.by value .by of the earlier sample.
#' @importFrom SummarizedExperiment colData
.change_table <- function(x, time.col, group, assay.type, detection,
        group.by, max.gap) {
    .check_longitudinal(x, time.col, group)
    .check_column(x, group.by)
    mat <- .get_assay(x, assay.type)
    pairs <- .consecutive_samples(x, time.col, group, max.gap)
    before <- mat[, pairs$from, drop = FALSE]
    after <- mat[, pairs$to, drop = FALSE]
    ok <- which(.present(before, detection) & .present(after, detection))
    if (length(ok) == 0L) {
        stop("No feature is above 'detection' in two consecutive samples.",
            call. = FALSE)
    }
    j <- col(before)[ok]
    by <- if (is.null(group.by)) "" else
        as.character(colData(x)[[group.by]])[pairs$from[j]]
    data.frame(.change = log(after[ok] / before[ok]),
        .interval = pairs$interval[j], .by = by)
}

# Density of the values v in the bins between breaks, at the centres of
# the bins with values.
.change_density <- function(v, breaks) {
    n <- tabulate(findInterval(v, breaks, all.inside = TRUE),
        length(breaks) - 1L)
    mids <- (breaks[-1L] + breaks[-length(breaks)]) / 2
    data.frame(.x = mids, .y = n / (length(v) * diff(breaks)))[n > 0L, ,
        drop = FALSE]
}

# Laplace density fitted by maximum likelihood over the range of v; .label
# gives the scale.
#' @importFrom stats median
.fit_laplace <- function(v) {
    m <- median(v)
    b <- mean(abs(v - m))
    if (b == 0) {
        return(NULL)
    }
    x <- seq(min(v), max(v), length.out = 200L)
    data.frame(.x = x, .y = exp(-abs(x - m) / b) / (2 * b),
        .label = sprintf("%.2f", b))
}

#' @importFrom ggplot2 ggplot aes .data geom_point geom_line labs
#'   scale_y_log10
.plot_changes <- function(df, bins, curves, time.col, group.by, ...) {
    interval <- unique(range(df$.interval))
    subtitle <- sprintf("%d changes, %s intervals %s", nrow(df), time.col,
        paste(format(interval), collapse = " to "))
    p <- ggplot(bins, aes(x = .data$.x, y = .data$.y,
            colour = if (!is.null(group.by)) .data$.by)) +
        geom_point(...)
    if (!is.null(curves)) {
        fits <- unique(curves[c(".by", ".label")])
        p <- p + geom_line(data = curves, linetype = "dashed")
        subtitle <- paste0(subtitle, "; Laplace scale ", paste(
            if (is.null(group.by)) fits$.label else
                paste(fits$.by, fits$.label), collapse = ", "))
    }
    p + scale_y_log10() +
        .scale_categorical("colour") +
        labs(x = "log(x[t+1] / x[t])", y = "density", colour = group.by,
            subtitle = subtitle) +
        .theme_tempoPlot()
}

# Candidate for utils.R (plotTurnover() pairs samples the same way).
# Consecutive samples of each subject, not across a gap larger than
# max.gap: column indices from and to of x and their time interval.
#' @importFrom SummarizedExperiment colData
.consecutive_samples <- function(x, time.col, group, max.gap) {
    cd <- as.data.frame(colData(x)[, c(time.col, group), drop = FALSE],
        optional = TRUE)
    cd$.i <- seq_len(ncol(x))
    cd <- .break_gaps(.order_time(cd, time.col, group), time.col, group,
        max.gap)
    subject <- if (is.null(group)) rep("", nrow(cd)) else
        as.character(cd[[group]])
    i <- .next_in_subject(subject, segment = cd$.segment)
    time <- .time_numeric(cd[[time.col]])
    list(from = cd$.i[i], to = cd$.i[i + 1L],
        interval = time[i + 1L] - time[i])
}
