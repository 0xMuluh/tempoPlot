#' Plot the noise colour of features per subject
#'
#' \code{plotNoise()} classifies the series of each feature in each subject
#' by the slope of its power spectrum as white, pink, brown or black noise,
#' and draws the share of features of each colour as a stacked bar per
#' subject.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date}
#'   or \code{POSIXct}.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of the
#'   features. (Default: \code{"counts"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_bar]{geom_bar}}.
#'
#' @details
#' Each series is the Lomb-Scargle periodogram of \code{\link{plotSpectrum}}
#' (linear trend removed), with a straight line fitted on log-log axes. As
#' in Faust et al. (2018, doi:10.1186/s40168-018-0496-2), a slope below
#' -2.25 is black noise, from -2.25 to -1.75 brown, from -1.75 to -0.5
#' pink, and above white. White noise means no dependence on earlier
#' samples; the steeper the slope, the longer the memory. Faust et al.
#' take the steepest slope of a smoothing spline, which is this straight
#' line for series of up to 100 samples. Faust et al. used relative
#' abundances of evenly sampled series; uneven sampling flattens the
#' slope towards white noise.
#'
#' A feature with fewer than four non-missing samples in a subject, or
#' with the same value in all of them, has no colour there and is left
#' out. The plot data hold the slope and colour of each feature.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @seealso
#' \code{\link{plotSpectrum}} with \code{fit = "power"} for the spectrum
#' and slope of one feature.
#'
#' @name plotNoise
#'
#' @examples
#' # Hourly samples with reads of four vessels of an artificial gut, hours
#' # 476 to 592, one sample per hour, and the 50 most abundant taxa
#' data("SilvermanAGutData", package = "miaTime")
#' tse <- SilvermanAGutData
#' tse <- tse[, tse$SampleType == "Hourly" &
#'     tse$Normal_Noise_Sample == "Normal" & tse$Time.hr %in% 476:592 &
#'     colSums(SummarizedExperiment::assay(tse)) > 0]
#' tse <- tse[, !duplicated(paste(tse$Vessel, tse$Time.hr))]
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' tse <- tse[mia::getTop(tse, top = 50), ]
#' plotNoise(tse, time.col = "Time.hr", group = "Vessel",
#'     assay.type = "relabundance")
NULL

#' @rdname plotNoise
#' @importFrom methods setMethod
#' @importFrom stats var
#' @export
setMethod("plotNoise", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts", ...) {
        design <- .spectral_design(x, time.col, group)
        mat <- .get_assay(x, assay.type)
        df <- do.call(rbind, lapply(rownames(mat), function(feature) {
            series <- .subject_series(design, mat[feature, ])
            series <- series[vapply(series, function(s) var(s$v) > 0,
                logical(1))]
            if (length(series) == 0L) {
                return(NULL)
            }
            data.frame(.subject = names(series), .feature = feature,
                slope = vapply(series, function(s)
                    .spectral_slope(.lomb_scargle(s$t, s$v)), numeric(1)))
        }))
        if (is.null(df)) {
            stop("No feature varies over four or more samples of a ",
                "subject.", call. = FALSE)
        }
        df$noise <- .noise_colour(df$slope)
        .plot_noise(df, group, ...)
    }
)

# Noise colour of spectral slopes, with the thresholds of Faust et al.
# (2018); levels from white to black.
.noise_colour <- function(slope) {
    colour <- cut(slope, c(-Inf, -2.25, -1.75, -0.5, Inf),
        labels = c("black", "brown", "pink", "white"), right = FALSE)
    factor(colour, levels = rev(levels(colour)))
}

#' @importFrom ggplot2 ggplot aes .data geom_bar labs scale_fill_manual
.plot_noise <- function(df, group, ...) {
    # Bars in the colour they name; the outline keeps white visible.
    ggplot(df, aes(x = .data$.subject, fill = .data$noise)) +
        # show.legend keeps the keys of colours that no feature has.
        geom_bar(position = "fill", colour = "grey30", show.legend = TRUE,
            ...) +
        scale_fill_manual(values = c(white = "white", pink = "pink",
            brown = "sienna4", black = "black"), drop = FALSE) +
        labs(x = group, y = "share of features", fill = "noise") +
        .theme_tempoPlot()
}
