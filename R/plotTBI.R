#' Plot losses against gains between two time points
#'
#' \code{plotTBI()} splits the change of each subject between two time
#' points into losses and gains, the components of the temporal
#' beta-diversity index (TBI), and draws one point per subject. Subjects
#' above the 1:1 line gained more than they lost.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param time.point \code{Vector} of two values of \code{time.col}: the
#'   first and the second time point. Subjects without a sample at both are
#'   left out.
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of
#'   abundances. (Default: \code{"counts"})
#'
#' @param method \code{Character scalar}. Splits the Sorensen
#'   dissimilarity of presence and absence (\code{"sorensen"}) or the
#'   percentage difference (Bray-Curtis) of abundances (\code{"bray"}).
#'   (Default: \code{"sorensen"})
#'
#' @param colour.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that colours the points, taken from the sample at
#'   the first time point. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' For abundances \eqn{y_1} and \eqn{y_2}, \eqn{A = \sum \min(y_1, y_2)},
#' losses \eqn{B = \sum \max(y_1 - y_2, 0)} and gains
#' \eqn{C = \sum \max(y_2 - y_1, 0)}. Both axes are divided by
#' \eqn{2A + B + C}, so that losses and gains add up to the dissimilarity.
#' With \code{"sorensen"}, the abundances are presence (1) and absence (0),
#' and losses and gains count the features lost and gained.
#'
#' \code{"bray"} needs absolute abundances, such as counts scaled by qPCR
#' or flow cytometry. \eqn{B - C} is the difference of the sample totals,
#' so on relative abundances losses always equal gains, and on sequencing
#' counts they differ only by library size.
#'
#' Legendre's permutation test of subjects that changed in exceptional
#' ways is not included.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Legendre P (2019). A temporal beta-diversity index to identify sites
#' that have changed in exceptional ways in space-time surveys.
#' \emph{Ecology and Evolution} 9(6):3500-3514.
#'
#' @name plotTBI
#'
#' @examples
#' # Genera lost and gained between the first and the second sample of
#' # each subject
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' plotTBI(tse, time.col = "visit", group = "subject",
#'     time.point = c("1", "2"), colour.by = "sex")
NULL

#' @rdname plotTBI
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotTBI", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, time.point, assay.type = "counts",
            method = c("sorensen", "bray"), colour.by = NULL, ...) {
        method <- match.arg(method)
        .require_column(x, group)
        .check_longitudinal(x, time.col, group)
        .check_column(x, colour.by)
        mat <- .get_assay(x, assay.type)
        if (method == "sorensen") {
            mat <- .present(mat) + 0
        }
        i <- .time_point_index(x, time.col, time.point)
        subject <- as.character(colData(x)[[group]])
        i1 <- i[[1L]]
        i2 <- i[[2L]]
        both <- intersect(subject[i1], subject[i2])
        if (length(both) == 0L) {
            stop("No subject has samples at both values of 'time.point'.",
                call. = FALSE)
        }
        i1 <- i1[match(both, subject[i1])]
        i2 <- i2[match(both, subject[i2])]
        y1 <- mat[, i1, drop = FALSE]
        y2 <- mat[, i2, drop = FALSE]
        a <- colSums(pmin(y1, y2), na.rm = TRUE)
        b <- colSums(pmax(y1 - y2, 0), na.rm = TRUE)
        g <- colSums(pmax(y2 - y1, 0), na.rm = TRUE)
        den <- 2 * a + b + g
        df <- data.frame(.subject = both,
            losses = ifelse(den > 0, b / den, 0),
            gains = ifelse(den > 0, g / den, 0))
        if (!is.null(colour.by)) {
            df[[colour.by]] <- colData(x)[[colour.by]][i1]
        }
        .plot_tbi(df, time.point, colour.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_point geom_abline labs
#'   coord_equal
.plot_tbi <- function(df, time.point, colour.by, ...) {
    lim <- c(0, max(c(df$losses, df$gains), 1e-12))
    ggplot(df, aes(x = .data$losses, y = .data$gains,
            colour = if (!is.null(colour.by)) .data[[colour.by]])) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
            colour = "grey40") +
        geom_point(...) +
        .scale_by(df, colour.by, "colour") +
        coord_equal(xlim = lim, ylim = lim) +
        labs(x = paste("losses from", time.point[1L], "to", time.point[2L]),
            y = paste("gains from", time.point[1L], "to", time.point[2L]),
            colour = colour.by) +
        .theme_tempoPlot()
}
