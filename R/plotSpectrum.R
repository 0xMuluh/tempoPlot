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
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date}
#'   or \code{POSIXct}.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to
#'   retrieve \code{y} from when \code{y} is a feature name.
#'   (Default: \code{"counts"})
#'
#' @param fit \code{Character scalar}. Draws no fit (\code{"none"}) or a
#'   straight line on log-log axes (\code{"power"}) for each subject; the
#'   subtitle gives its slope. (Default: \code{"none"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' Frequencies are in cycles per \code{time.col} unit (per day for
#' \code{Date}, per hour for \code{POSIXct} times), from one cycle over
#' the series to half the median sampling rate. A linear trend is removed
#' first, and the first and last tenth of the series are tapered with a
#' cosine bell, as \code{\link[stats:spectrum]{spectrum}} does by default;
#' the taper reduces the leakage of power from low to high frequencies.
#' The power is normalised by the variance of the series. The
#' periodogram allows uneven sampling; for even sampling it equals the
#' classical periodogram. Subjects with fewer than four samples are left
#' out.
#'
#' The slope of the log-log fit gives the noise colour of the series:
#' about 0 for white noise, -1 for pink, -2 for brown and below for black
#' noise (Faust et al. 2018, doi:10.1186/s40168-018-0496-2).
#' \code{\link{plotNoise}} classifies the features of each subject by
#' this slope.
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
#'
#' # Hourly samples with reads of four vessels of an artificial gut, hours
#' # 476 to 592, one sample per hour: the slope of the spectrum of the most
#' # abundant taxon on log-log axes
#' data("SilvermanAGutData", package = "miaTime")
#' tse <- SilvermanAGutData
#' tse <- tse[, tse$SampleType == "Hourly" &
#'     tse$Normal_Noise_Sample == "Normal" & tse$Time.hr %in% 476:592 &
#'     colSums(SummarizedExperiment::assay(tse)) > 0]
#' tse <- tse[, !duplicated(paste(tse$Vessel, tse$Time.hr))]
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' plotSpectrum(tse, y = mia::getTop(tse, top = 1), time.col = "Time.hr",
#'     group = "Vessel", assay.type = "relabundance", fit = "power")
NULL

#' @rdname plotSpectrum
#' @importFrom methods setMethod
#' @export
setMethod("plotSpectrum", signature = c(x = "SummarizedExperiment"),
    function(x, y, time.col, group = NULL, assay.type = "counts",
            fit = c("none", "power"), ...) {
        fit <- match.arg(fit)
        series <- .subject_series(.spectral_design(x, time.col, group),
            unname(.resolve_y(x, y, assay.type)))
        if (length(series) == 0L) {
            stop("No subject has four or more samples of 'y'.",
                call. = FALSE)
        }
        df <- do.call(rbind, Map(function(s, id)
            cbind(.subject = id, .lomb_scargle(s$t, s$v)),
            series, names(series)))
        .plot_spectrum(df, y, time.col, group, fit, ...)
    }
)

# Times as numbers and subjects of the samples of x; time.col must be
# continuous.
#' @importFrom SummarizedExperiment colData
.spectral_design <- function(x, time.col, group) {
    if (.check_longitudinal(x, time.col, group) == "visit") {
        stop("'time.col' must be numeric, Date or POSIXct.",
            call. = FALSE)
    }
    list(t = .time_numeric(colData(x)[[time.col]]),
        subject = if (is.null(group)) rep("", ncol(x)) else
            as.character(colData(x)[[group]]))
}

# Series of value (one per sample) for each subject of design, ordered by
# time, without missing values. Subjects with fewer than four samples are
# left out.
.subject_series <- function(design, value) {
    keep <- !is.na(value)
    series <- split(data.frame(t = design$t[keep], v = value[keep]),
        design$subject[keep])
    series <- lapply(series, function(s) s[order(s$t), , drop = FALSE])
    series[vapply(series, nrow, integer(1)) >= 4L]
}

# Normalised Lomb-Scargle periodogram on a grid from one cycle over the
# series to half the median sampling rate.
#' @importFrom stats median var lm residuals
.lomb_scargle <- function(t, v, n = 200L) {
    span <- max(t) - min(t)
    # Remove a linear trend and taper the ends: a split cosine bell over
    # the first and last tenth of the span, as in stats::spec.taper().
    u <- (t - min(t)) / span
    taper <- (1 - cos(pi * pmin(1, 10 * u, 10 * (1 - u)))) / 2
    v <- residuals(lm(v ~ t)) * taper
    s2 <- var(v)
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

# Slope of a straight line through a periodogram on log-log axes.
#' @importFrom stats lm coef
.spectral_slope <- function(spectrum) {
    coef(lm(log10(power) ~ log10(frequency), data = spectrum))[[2L]]
}

#' @importFrom ggplot2 ggplot aes .data geom_line geom_smooth labs
#'   scale_x_log10 scale_y_log10
.plot_spectrum <- function(df, y, time.col, group, fit = "none", ...) {
    p <- ggplot(df, aes(x = .data$frequency, y = .data$power,
            colour = if (!is.null(group)) .data$.subject,
            group = .data$.subject)) +
        geom_line(...) +
        .scale_categorical("colour") +
        labs(x = paste("frequency (cycles per", time.col, "unit)"),
            y = paste("power of", y), colour = group) +
        .theme_tempoPlot()
    if (fit == "power") {
        # geom_smooth() fits on the log scales: the line of
        # .spectral_slope().
        slope <- vapply(split(df, df$.subject), .spectral_slope,
            numeric(1))
        label <- sprintf("slope %.2f", slope)
        p <- p + do.call(geom_smooth, c(list(method = "lm",
                formula = y ~ x, se = FALSE, linetype = "dashed"),
                if (is.null(group)) list(colour = "black"))) +
            scale_x_log10() +
            scale_y_log10() +
            labs(subtitle = paste(if (is.null(group)) label else
                paste0(names(slope), ": ", label), collapse = "; "))
    }
    p
}
