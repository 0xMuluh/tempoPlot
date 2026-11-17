#' Plot the accumulation of features over time
#'
#' \code{plotAccumulation()} counts the distinct features observed in each
#' subject up to each time point and draws one line per subject: a temporal
#' species accumulation curve. Curves that keep rising show subjects that
#' keep gaining new taxa; flat curves show a closed community.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of
#'   abundances. (Default: \code{"counts"})
#'
#' @param detection \code{Numeric scalar}. A feature is observed in a
#'   sample when its value is above this. (Default: \code{0})
#'
#' @param rank \code{Character scalar}. Specifies the column of
#'   \code{rowData(x)}, such as a taxonomic rank, whose levels are counted
#'   instead of features. (Default: \code{NULL})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that colours subjects and adds the mean curve of each
#'   group at each time; best with times shared by the subjects.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' Missing values count as not observed. A feature observed once stays
#' counted.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotAccumulation
#'
#' @examples
#' # Daily samples of four artificial guts
#' data("SilvermanAGutData", package = "miaTime")
#' tse <- SilvermanAGutData[, SilvermanAGutData$SampleType == "Daily"]
#' # Some vessels have two daily samples on a day; keep the first
#' tse <- tse[, !duplicated(paste(tse$Vessel, tse$DAY_ORDER))]
#' plotAccumulation(tse, time.col = "DAY_ORDER", group = "Vessel")
#'
#' # Genera instead of sequence variants
#' plotAccumulation(tse, time.col = "DAY_ORDER", group = "Vessel",
#'     rank = "Genus")
NULL

#' @rdname plotAccumulation
#' @importFrom SummarizedExperiment colData assay assayNames rowData
#' @importFrom stats ave
#' @importFrom methods setMethod
#' @export
setMethod("plotAccumulation", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            detection = 0, rank = NULL, group.by = NULL, ...) {
        .check_longitudinal(x, time.col, group)
        .time_type(x, time.col)
        .check_column(x, group.by)
        .check_column(x, rank, where = "rowData")
        if (!.is_one_name(assay.type) || !assay.type %in% assayNames(x)) {
            stop("'assay.type' must be one of assayNames(x).",
                call. = FALSE)
        }
        if (!(is.numeric(detection) && length(detection) == 1L)) {
            stop("'detection' must be a single number.", call. = FALSE)
        }
        present <- as.matrix(assay(x, assay.type)) > detection
        present[is.na(present)] <- FALSE
        if (!is.null(rank)) {
            levels <- as.character(rowData(x)[[rank]])
            keep <- !is.na(levels)
            present <- rowsum(present[keep, , drop = FALSE] + 0,
                levels[keep]) > 0
        }

        cols <- unique(c(time.col, group, group.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df$.i <- seq_len(ncol(x))
        df <- .order_time(df, time.col, group)
        subject <- if (is.null(group)) rep("", nrow(df)) else
            as.character(df[[group]])
        df$.subject <- subject
        df$.n <- unlist(lapply(split(df$.i, factor(subject,
                levels = unique(subject))), function(i) {
            # A feature seen in a sample stays seen in the later ones.
            seen <- present[, i, drop = FALSE]
            for (k in seq_len(ncol(seen))[-1L]) {
                seen[, k] <- seen[, k] | seen[, k - 1L]
            }
            colSums(seen)
        }), use.names = FALSE)
        .plot_accumulation(df, time.col, rank, group.by, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point stat_summary
#'   labs scale_x_discrete
.plot_accumulation <- function(df, time.col, rank, group.by, ...) {
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data$.n,
            colour = if (!is.null(group.by)) .data[[group.by]]))
    if (is.null(group.by)) {
        p <- p + geom_line(aes(group = .data$.subject), ...) +
            geom_point(size = 1)
    } else {
        p <- p + geom_line(aes(group = .data$.subject), alpha = 0.3, ...) +
            stat_summary(aes(group = .data[[group.by]]), fun = mean,
                geom = "line", linewidth = 1)
    }
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p + labs(x = time.col, colour = group.by,
        y = paste("cumulative", if (is.null(rank)) "features" else rank))
}
