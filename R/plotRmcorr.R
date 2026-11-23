#' Plot the repeated-measures correlation of two variables
#'
#' \code{plotRmcorr()} draws two variables of the same samples against each
#' other, coloured by subject, with one fitted line per subject. The lines
#' share the common within-subject slope, and the subtitle gives the
#' repeated-measures correlation.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param x.var,y.var \code{Character scalar}. Variables on the x and y
#'   axes: a column of \code{colData(x)}, a feature of \code{assay.type}, or
#'   \code{"<dimred>:<component>"}.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param assay.type \code{Character scalar}. Specifies the assay that holds
#'   features named in \code{x.var} or \code{y.var}.
#'   (Default: \code{"counts"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' The repeated-measures correlation of Bakdash and Marusich (2017,
#' doi:10.3389/fpsyg.2017.00456) is the association of the two variables
#' within subjects. It comes from the analysis of covariance
#' \code{lm(y.var ~ group + x.var)}:
#' \eqn{r_{rm} = \pm\sqrt{SS_x / (SS_x + SS_{error})}}, with the sign of
#' the common slope, and its p-value from the F test of \code{x.var} with
#' \eqn{L - N - 1} error degrees of freedom, for \eqn{L} samples of
#' \eqn{N} subjects. This equals the Pearson correlation of the two
#' variables after centring each on its subject means. Each line passes
#' through its subject's means and spans the subject's values of
#' \code{x.var}. Samples with a missing value are dropped.
#'
#' The test treats the samples of a subject as independent given the
#' subject, so serial correlation in densely sampled time series makes the
#' p-value too small. Zhou et al. (2019, doi:10.1038/s41586-019-1236-x,
#' Fig. 5b) use the plot to relate host and microbial measures within
#' participants. The colour legend is shown for up to six subjects; with
#' more, the colours cannot identify subjects and the legend is dropped.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Bakdash JZ, Marusich LR (2017). Repeated measures correlation.
#' \emph{Frontiers in Psychology} 8:456.
#'
#' @name plotRmcorr
#'
#' @examples
#' # Diversity against the relative abundance of a genus within subjects
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' plotRmcorr(tse, x.var = "Bacteroides vulgatus et rel.",
#'     y.var = "diversity", group = "subject", assay.type = "relabundance")
#'
#' # Two species within three bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "clr", pseudocount = 1)
#' plotRmcorr(tse, x.var = rownames(tse)[1L], y.var = rownames(tse)[2L],
#'     group = "StudyIdentifier", assay.type = "clr")
NULL

#' @rdname plotRmcorr
#' @importFrom SummarizedExperiment colData
#' @importFrom stats ave complete.cases
#' @importFrom methods setMethod
#' @export
setMethod("plotRmcorr", signature = c(x = "SummarizedExperiment"),
    function(x, x.var, y.var, group, assay.type = "counts", ...) {
        .require_column(x, group)
        df <- data.frame(
            .subject = colData(x)[[group]],
            .x = unname(.resolve_y(x, x.var, assay.type, "x.var")),
            .y = unname(.resolve_y(x, y.var, assay.type, "y.var")))
        df <- df[complete.cases(df), , drop = FALSE]
        df$.subject <- factor(df$.subject)
        fit <- .rmcorr(df$.x, df$.y, df$.subject)
        df$.fit <- ave(df$.y, df$.subject) +
            fit$slope * (df$.x - ave(df$.x, df$.subject))
        .plot_rmcorr(df, fit, x.var, y.var, group, ...)
    }
)

# Repeated-measures correlation r of x and y within subjects (a factor
# without unused levels): the Pearson correlation after centring both on
# the subject means, which is the partial correlation of y ~ subject + x.
# Also the common slope and the t test of r with length(x) - subjects - 1
# degrees of freedom, which equals the F test of x in that model.
#' @importFrom stats ave pt
.rmcorr <- function(x, y, subject) {
    xc <- x - ave(x, subject)
    yc <- y - ave(y, subject)
    df <- length(x) - nlevels(subject) - 1L
    if (df < 1L || sum(xc^2) == 0 || sum(yc^2) == 0) {
        stop("'x.var' and 'y.var' must have values that vary within ",
            "subjects of 'group'.", call. = FALSE)
    }
    r <- sum(xc * yc) / sqrt(sum(xc^2) * sum(yc^2))
    t <- r * sqrt(df / (1 - r^2))
    list(r = r, df = df, p = 2 * pt(-abs(t), df),
        slope = sum(xc * yc) / sum(xc^2))
}

#' @importFrom ggplot2 ggplot aes .data geom_point geom_line labs
.plot_rmcorr <- function(df, fit, x.var, y.var, group, ...) {
    legend <- if (nlevels(df$.subject) > length(.okabe_ito)) "none" else
        "legend"
    ggplot(df, aes(x = .data$.x, colour = .data$.subject)) +
        geom_line(aes(y = .data$.fit, group = .data$.subject)) +
        geom_point(aes(y = .data$.y), ...) +
        .scale_categorical("colour", guide = legend) +
        labs(x = x.var, y = y.var, colour = group,
            subtitle = sprintf("r_rm = %.2f, df = %d, p = %.2g", fit$r,
                fit$df, fit$p)) +
        .theme_tempoPlot()
}
