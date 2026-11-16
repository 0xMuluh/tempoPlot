#' Plot whether subjects are identified by their own earlier sample
#'
#' \code{plotFingerprint()} compares, for each subject, the sample at a
#' query time with all samples at a reference time. A subject is
#' identified when its own reference sample is closer than the reference
#' sample of any other subject.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param dissimilarity \code{Character scalar}. Specifies the
#'   sample by sample dissimilarity in \code{metadata(x)}, as stored by
#'   \code{\link[mia:getDissimilarity]{mia::addDissimilarity()}}.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param reference.time,query.time Values of \code{time.col}: the samples
#'   to identify from (\code{reference.time}) and the samples to identify
#'   (\code{query.time}).
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Each point is a subject with samples at both times: the dissimilarity
#' of its query sample to its own reference sample (x) and to the nearest
#' reference sample of another subject (y). Points above the 1:1 line are
#' identified. The subtitle gives the share of subjects identified.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Franzosa EA, et al. (2015). Identifying personal microbiomes using
#' metagenomic codes. \emph{PNAS} 112(22):E2930-E2938.
#'
#' @name plotFingerprint
#'
#' @examples
#' # Is each subject's second sample closest to its own first sample?
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' tse <- tse[, tse$visit %in% c("1", "2")]
#' tse <- mia::addDissimilarity(tse, method = "bray", assay.type = "counts")
#' plotFingerprint(tse, dissimilarity = "bray", time.col = "visit",
#'     group = "subject", reference.time = "1", query.time = "2")
NULL

#' @rdname plotFingerprint
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotFingerprint", signature = c(x = "SummarizedExperiment"),
    function(x, dissimilarity, time.col, group, reference.time, query.time,
            ...) {
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        .time_type(x, time.col)
        d <- .get_dissimilarity(x, dissimilarity)
        time <- as.character(colData(x)[[time.col]])
        for (arg in c("reference.time", "query.time")) {
            v <- get(arg)
            if (length(v) != 1L || !as.character(v) %in% time) {
                stop("'", arg, "' must be one value of 'time.col'.",
                    call. = FALSE)
            }
        }
        subject <- as.character(colData(x)[[group]])
        ref <- which(time == as.character(reference.time))
        qry <- which(time == as.character(query.time))
        qry <- qry[subject[qry] %in% subject[ref]]
        if (length(qry) == 0L || length(unique(subject[ref])) < 2L) {
            stop("Subjects need samples at both times, and at least two ",
                "subjects at 'reference.time'.", call. = FALSE)
        }
        df <- do.call(rbind, lapply(qry, function(q) {
            own <- ref[subject[ref] == subject[q]]
            other <- ref[subject[ref] != subject[q]]
            data.frame(.subject = subject[q], own = d[q, own],
                other = min(d[q, other]))
        }))
        df$identified <- df$own < df$other
        .plot_fingerprint(df, dissimilarity, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_point geom_abline labs
#'   scale_colour_manual coord_equal
.plot_fingerprint <- function(df, dissimilarity, ...) {
    lim <- range(c(0, df$own, df$other))
    ggplot(df, aes(x = .data$own, y = .data$other,
            colour = .data$identified)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
            colour = "grey40") +
        geom_point(...) +
        scale_colour_manual(values = c(`TRUE` = "#2c6fbb",
            `FALSE` = "#b2182b"), labels = c(`TRUE` = "yes",
            `FALSE` = "no")) +
        coord_equal(xlim = lim, ylim = lim) +
        labs(x = paste(dissimilarity, "to own reference sample"),
            y = paste(dissimilarity, "to nearest other subject"),
            colour = "identified",
            subtitle = sprintf("%d of %d subjects identified (%.0f%%)",
                sum(df$identified), nrow(df), 100 * mean(df$identified)))
}
