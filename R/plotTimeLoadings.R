#' Plot temporal loadings of a tensor decomposition
#'
#' \code{plotTimeLoadings()} draws the time loading function of each
#' component of a tensor decomposition stored in \code{metadata(x)}, such
#' as TEMPTED or PARAFAC: how strongly the component is expressed at each
#' time.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param source \code{Character scalar}. Specifies the decomposition in
#'   \code{metadata(x)}: a result of \code{tempted::tempted_all()} or
#'   \code{tempted::tempted()}, a PARAFAC model with the loadings in
#'   \code{Fac} (as from \code{parafac4microbiome::parafac()}), whose third
#'   mode is time, or a matrix of time points by components with the times
#'   as row names or in a column \code{time}.
#'
#' @param ncomponents \code{Integer scalar}. Number of components to draw,
#'   from the first. If \code{NULL}, draws all. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' For TEMPTED, the loadings are functions evaluated on a time grid, and
#' the legend gives the share of variance of each component
#' (\code{r_square}). For PARAFAC, the time mode has one loading per time
#' point; set its row names to the times, otherwise points are numbered in
#' their order. The sign of a component is arbitrary: a component and its
#' negative describe the same pattern. Feature and subject loadings are
#' drawn by \code{\link[miaViz:plotLoadings]{miaViz::plotLoadings()}} and
#' \code{\link{plotSubjectScatter}}.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Shi P, et al. (2024). TEMPTED: time-informed dimensionality reduction
#' for longitudinal microbiome studies. \emph{Genome Biology} 25:317.
#'
#' @name plotTimeLoadings
#'
#' @examples
#' # TEMPTED of the species in three bioreactors
#' if (requireNamespace("tempted", quietly = TRUE)) {
#'     data("minimalgut", package = "miaTime")
#'     tse <- minimalgut
#'     counts <- t(SummarizedExperiment::assay(tse, "counts"))
#'     S4Vectors::metadata(tse)$tempted <- tempted::tempted_all(counts,
#'         timepoint = tse$Time_hr_num,
#'         subjectID = as.character(tse$StudyIdentifier), r = 2,
#'         threshold = 1, pseudo = 0.5, do_ratio = FALSE)
#'     plotTimeLoadings(tse, source = "tempted")
#' }
NULL

#' @rdname plotTimeLoadings
#' @importFrom S4Vectors metadata
#' @importFrom methods setMethod
#' @export
setMethod("plotTimeLoadings", signature = c(x = "SummarizedExperiment"),
    function(x, source, ncomponents = NULL, ...) {
        if (!.is_one_name(source) || !source %in% names(metadata(x))) {
            stop("'source' must name a decomposition in metadata(x).",
                call. = FALSE)
        }
        res <- .time_loadings(metadata(x)[[source]])
        k <- ncol(res$loadings)
        if (!is.null(ncomponents)) {
            if (!(is.numeric(ncomponents) && length(ncomponents) == 1L &&
                    ncomponents >= 1 && ncomponents <= k)) {
                stop("'ncomponents' must be between 1 and ", k, ".",
                    call. = FALSE)
            }
            k <- as.integer(ncomponents)
        }
        labels <- res$labels[seq_len(k)]
        df <- data.frame(
            .time = rep(res$time, times = k),
            .component = factor(rep(labels, each = length(res$time)),
                levels = labels),
            .value = as.vector(res$loadings[, seq_len(k), drop = FALSE]))
        .plot_time_loadings(df, res$discrete, ...)
    }
)

# Time loadings, times and component labels of a decomposition.
.time_loadings <- function(obj) {
    if (is.list(obj) && all(c("Phi_hat", "time_Phi") %in% names(obj))) {
        mat <- as.matrix(obj$Phi_hat)
        labels <- colnames(mat) %||% paste0("PC", seq_len(ncol(mat)))
        if (length(obj$r_square) == ncol(mat)) {
            labels <- sprintf("%s (%.0f%%)", labels, 100 * obj$r_square)
        }
        return(list(loadings = mat, time = as.numeric(obj$time_Phi),
            labels = labels, discrete = FALSE))
    }
    if (is.list(obj) && !is.data.frame(obj) && "Fac" %in% names(obj) &&
            length(obj$Fac) >= 3L) {
        mat <- as.matrix(obj$Fac[[3L]])
        time <- suppressWarnings(as.numeric(rownames(mat)))
        if (is.null(rownames(mat)) || anyNA(time)) {
            time <- seq_len(nrow(mat))
        }
        return(list(loadings = mat, time = time,
            labels = paste0("Component ", seq_len(ncol(mat))),
            discrete = TRUE))
    }
    if (is.matrix(obj) || is.data.frame(obj)) {
        df <- as.data.frame(obj)
        if ("time" %in% names(df)) {
            time <- df$time
            df$time <- NULL
        } else {
            time <- suppressWarnings(as.numeric(rownames(df)))
        }
        if (anyNA(time) || !all(vapply(df, is.numeric, logical(1)))) {
            stop("A loading matrix needs numeric times, as row names or ",
                "in a column 'time', and numeric loadings.", call. = FALSE)
        }
        return(list(loadings = as.matrix(df), time = as.numeric(time),
            labels = names(df) %||% paste0("Component ", seq_len(ncol(df))),
            discrete = FALSE))
    }
    stop("'source' is not a TEMPTED result, a PARAFAC model or a loading ",
        "matrix.", call. = FALSE)
}

#' @importFrom ggplot2 ggplot aes .data geom_line geom_point geom_hline
#'   labs
.plot_time_loadings <- function(df, discrete, ...) {
    p <- ggplot(df, aes(x = .data$.time, y = .data$.value,
            colour = .data$.component)) +
        geom_hline(yintercept = 0, colour = "grey60") +
        geom_line(...)
    if (discrete) {
        p <- p + geom_point(size = 1)
    }
    p + labs(x = "time", y = "loading", colour = "component")
}
