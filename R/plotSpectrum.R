#' Plot the power spectrum of a series
#'
#' \code{plotSpectrum()} draws the Lomb-Scargle periodogram of a variable,
#' one curve per subject. Peaks show periodic changes, such as daily
#' rhythms.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param y \code{Character scalar}. Specifies the variable: a column of
#'   \code{colData(x)}, a row of \code{assay(x, assay.type)}, or
#'   \code{"<dimred>:<component>"}.
#'
#' @param time.col \code{Character scalar}. Specifies the numeric or
#'   \code{Date} column from \code{colData(x)} that contains sampling times.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to
#'   retrieve \code{y} from when \code{y} is a feature name.
#'   (Default: \code{"counts"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' Frequencies are in cycles per \code{time.col} unit (per day for
#' \code{Date} times), from one cycle over the series to half the median
#' sampling rate. A linear trend is removed first. The power is normalised
#' by the variance of the series. The periodogram allows uneven sampling;
#' for even sampling it equals the classical periodogram. Subjects with
#' fewer than four samples are left out.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Lomb NR (1976). Least-squares frequency analysis of unequally spaced
#' data. \emph{Astrophysics and Space Science} 39:447-462.
#'
#' Scargle JD (1982). Studies in astronomical time series analysis. II.
#' \emph{The Astrophysical Journal} 263:835-853.
#'
#' @name plotSpectrum
#'
#' @examples
#' # Acetate in three bioreactors fasted once a day: peaks near 1/24 per
#' # hour; the lowest frequencies carry the slow rise and fall of acetate
#' data("minimalgut", package = "miaTime")
#' plotSpectrum(minimalgut, y = "Acetate", time.col = "Time_hr_num",
#'     group = "StudyIdentifier")
NULL

#' @rdname plotSpectrum
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median var lm residuals
#' @importFrom methods setMethod
#' @export
setMethod("plotSpectrum", signature = c(x = "SummarizedExperiment"),
    function(x, y, time.col, group = NULL, assay.type = "counts",
            ...) {
        if (.check_longitudinal(x, time.col, group) == "visit") {
            stop("'time.col' must be numeric or Date.", call. = FALSE)
        }
        value <- unname(.resolve_y(x, y, assay.type))
        time <- as.numeric(colData(x)[[time.col]])
        subject <- if (is.null(group)) rep("", ncol(x)) else
            as.character(colData(x)[[group]])
        keep <- !is.na(value)
        series <- split(data.frame(t = time[keep], v = value[keep]),
            subject[keep])
        series <- series[vapply(series, nrow, integer(1)) >= 4L]
        if (length(series) == 0L) {
            stop("No subject has four or more samples of 'y'.",
                call. = FALSE)
        }
        df <- do.call(rbind, lapply(seq_along(series), function(i) {
            s <- series[[i]]
            s <- s[order(s$t), , drop = FALSE]
            cbind(.subject = names(series)[i], .lomb_scargle(s$t, s$v))
        }))
        .plot_spectrum(df, y, time.col, group, ...)
    }
)

# Normalised Lomb-Scargle periodogram on a grid from one cycle over the
# series to half the median sampling rate.
.lomb_scargle <- function(t, v, n = 200L) {
    # Remove a linear trend.
    v <- residuals(lm(v ~ t))
    s2 <- var(v)
    span <- max(t) - min(t)
    fmax <- 0.5 / median(diff(t))
    freq <- seq(1 / span, fmax, length.out = n)
    power <- vapply(freq, function(f) {
        w <- 2 * pi * f
        tau <- atan2(sum(sin(2 * w * t)), sum(cos(2 * w * t))) / (2 * w)
        c <- cos(w * (t - tau))
        s <- sin(w * (t - tau))
        (sum(v * c)^2 / sum(c^2) + sum(v * s)^2 / sum(s^2)) / (2 * s2)
    }, numeric(1))
    data.frame(frequency = freq, power = power)
}

#' @importFrom ggplot2 ggplot aes .data geom_line labs
.plot_spectrum <- function(df, y, time.col, group, ...) {
    ggplot(df, aes(x = .data$frequency, y = .data$power,
            colour = if (!is.null(group)) .data$.subject,
            group = .data$.subject)) +
        geom_line(...) +
        .scale_categorical("colour") +
        labs(x = paste("frequency (cycles per", time.col, "unit)"),
            y = paste("power of", y), colour = group) +
        .theme_tempoPlot()
}
