#' Plot dissimilarities between the samples of each subject
#'
#' \code{plotRecurrence()} draws, for each subject, a time by time matrix
#' of the dissimilarities between its samples. With \code{threshold}, it
#' draws a recurrence plot: the pairs of samples closer than the threshold.
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
#' @param subjects \code{Character vector}. Subjects to draw, one panel
#'   each. If \code{NULL}, draws all. (Default: \code{NULL})
#'
#' @param threshold \code{Numeric scalar}. Marks the pairs of samples with
#'   a dissimilarity of at most this. If \code{NULL}, colours by
#'   dissimilarity. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_rect]{geom_rect}}.
#'
#' @details
#' Each sample covers the time from halfway to its previous sample to
#' halfway to its next one, so the axes keep the spacing of the sampling
#' times. Visits are evenly spaced.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotRecurrence
#'
#' @examples
#' # Three bioreactors under the same perturbations
#' data("minimalgut", package = "miaTime")
#' tse <- mia::addDissimilarity(minimalgut, method = "bray",
#'     assay.type = "counts")
#' plotRecurrence(tse, dissimilarity = "bray", time.col = "Time.hr",
#'     group = "StudyIdentifier")
#'
#' # Pairs of samples within a Bray-Curtis dissimilarity of 0.2
#' plotRecurrence(tse, dissimilarity = "bray", time.col = "Time.hr",
#'     group = "StudyIdentifier", subjects = "Bioreactor A",
#'     threshold = 0.2)
NULL

#' @rdname plotRecurrence
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotRecurrence", signature = c(x = "SummarizedExperiment"),
    function(x, dissimilarity, time.col, group, subjects = NULL,
            threshold = NULL, ...) {
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        type <- .time_type(x, time.col)
        d <- .get_dissimilarity(x, dissimilarity)
        if (!is.null(threshold) && !(is.numeric(threshold) &&
                length(threshold) == 1L && is.finite(threshold))) {
            stop("'threshold' must be NULL or a single number.",
                call. = FALSE)
        }
        subject <- as.character(colData(x)[[group]])
        if (is.null(subjects)) {
            subjects <- unique(subject)
        } else if (!is.character(subjects) ||
                !all(subjects %in% subject)) {
            stop("'subjects' must be values of the 'group' column.",
                call. = FALSE)
        }

        time <- colData(x)[[time.col]]
        tnum <- .time_numeric(time)
        df <- do.call(rbind, lapply(subjects, function(s) {
            i <- which(subject == s)
            i <- i[order(tnum[i])]
            edges <- if (type == "visit") {
                list(lower = tnum[i] - 0.5, upper = tnum[i] + 0.5)
            } else {
                .tile_edges(tnum[i])
            }
            k <- expand.grid(a = seq_along(i), b = seq_along(i))
            data.frame(.subject = s,
                xmin = edges$lower[k$a], xmax = edges$upper[k$a],
                ymin = edges$lower[k$b], ymax = edges$upper[k$b],
                dissimilarity = d[cbind(i[k$a], i[k$b])])
        }))
        df$.subject <- factor(df$.subject, levels = subjects)
        names(df)[1L] <- group
        if (type == "date") {
            for (col in c("xmin", "xmax", "ymin", "ymax")) {
                df[[col]] <- as.Date(df[[col]], origin = "1970-01-01")
            }
        }
        levels <- if (type == "visit") levels(time)
        .plot_recurrence(df, dissimilarity, time.col, group, threshold,
            levels, ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_rect labs facet_wrap vars
#'   scale_fill_gradient scale_fill_manual scale_x_continuous
#'   scale_y_continuous
.plot_recurrence <- function(df, dissimilarity, time.col, group, threshold,
        levels, ...) {
    p <- ggplot(df, aes(xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax))
    if (is.null(threshold)) {
        p <- p + geom_rect(aes(fill = .data$dissimilarity), ...) +
            scale_fill_gradient(low = "#2c6fbb", high = "grey95") +
            labs(fill = dissimilarity)
    } else {
        df$recurrent <- df$dissimilarity <= threshold
        p <- p + geom_rect(aes(fill = .data$recurrent), data = df, ...) +
            scale_fill_manual(values = c(`TRUE` = "#2c6fbb",
                `FALSE` = "grey95"), labels = c(`TRUE` = "yes",
                `FALSE` = "no")) +
            labs(fill = paste(dissimilarity, "<=", threshold))
    }
    if (!is.null(levels)) {
        breaks <- seq_along(levels)
        p <- p + scale_x_continuous(breaks = breaks, labels = levels) +
            scale_y_continuous(breaks = breaks, labels = levels)
    }
    p + facet_wrap(vars(.data[[group]]), scales = "free") +
        labs(x = time.col, y = time.col)
}
