#' Plot two subject-level statistics against each other
#'
#' \code{plotSubjectScatter()} draws one point per subject from two
#' statistics that are constant within each subject, such as mean diversity
#' and its variability or resistance and resilience, or from two components
#' of the subject loadings of a tensor decomposition such as TEMPTED.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param x.stat,y.stat \code{Character scalar}. Specify the numeric columns
#'   from \code{colData(x)} drawn on the x and y axes. Not used with
#'   \code{source}. (Default: \code{NULL})
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param colour.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)}, constant within subjects, that colours the points.
#'   (Default: \code{NULL})
#'
#' @param quadrants \code{Logical scalar} or \code{Numeric vector} of
#'   length 2. Draws reference lines that split the plot into quadrants: at
#'   the medians of the two statistics (\code{TRUE}) or at the given x and y
#'   values. (Default: \code{FALSE})
#'
#' @param source \code{Character scalar}. Specifies a decomposition in
#'   \code{metadata(x)} whose subject loadings are drawn instead of
#'   \code{x.stat} and \code{y.stat}: a result of
#'   \code{tempted::tempted_all()} (\code{A_hat}), a PARAFAC model whose
#'   first mode is subjects (\code{Fac}, with subjects as row names), or a
#'   matrix of subjects by components with subjects as row names.
#'   (Default: \code{NULL})
#'
#' @param ncomponents \code{Numeric scalar} or \code{vector}. The two
#'   components of \code{source} on the x and y axes: \code{2} for the
#'   first two, or two components by index or name. (Default: \code{2})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#'
#' @details
#' Subjects with a missing statistic are dropped. Subject loadings are
#' matched to \code{group} by subject name; \code{colour.by} must be
#' constant within subjects.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotSubjectScatter
#'
#' @examples
#' # Mean diversity against its standard deviation over each subject's
#' # samples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
#' tse$mean <- ave(tse$diversity, tse$subject, FUN = mean)
#' tse$sd <- ave(tse$diversity, tse$subject, FUN = sd)
#' plotSubjectScatter(tse, x.stat = "mean", y.stat = "sd",
#'     group = "subject", colour.by = "sex", quadrants = TRUE)
#'
#' # TEMPTED subject loadings of the bioreactors
#' if (requireNamespace("tempted", quietly = TRUE)) {
#'     data("minimalgut", package = "miaTime")
#'     tse <- minimalgut
#'     counts <- t(SummarizedExperiment::assay(tse, "counts"))
#'     S4Vectors::metadata(tse)$tempted <- tempted::tempted_all(counts,
#'         timepoint = tse$Time_hr_num,
#'         subjectID = as.character(tse$StudyIdentifier), r = 2,
#'         threshold = 1, pseudo = 0.5, do_ratio = FALSE)
#'     plotSubjectScatter(tse, group = "StudyIdentifier",
#'         colour.by = "StudyIdentifier", source = "tempted")
#' }
NULL

#' @rdname plotSubjectScatter
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom stats median
#' @importFrom methods setMethod
#' @export
setMethod("plotSubjectScatter", signature = c(x = "SummarizedExperiment"),
    function(x, x.stat = NULL, y.stat = NULL, group, colour.by = NULL,
            quadrants = FALSE, source = NULL, ncomponents = 2, ...) {
        if (is.null(group)) {
            stop("'group' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_column(x, group)
        .check_column(x, colour.by)
        if (!(isTRUE(quadrants) || isFALSE(quadrants) ||
                (is.numeric(quadrants) && length(quadrants) == 2L))) {
            stop("'quadrants' must be TRUE, FALSE or two numbers.",
                call. = FALSE)
        }
        df <- if (is.null(source)) {
            .subject_stats(x, x.stat, y.stat, group, colour.by)
        } else {
            .subject_loading_table(x, source, ncomponents, group, colour.by)
        }
        labels <- attr(df, "labels")
        df <- df[!is.na(df$.x) & !is.na(df$.y), , drop = FALSE]
        if (nrow(df) == 0L) {
            stop("No subject has both statistics.", call. = FALSE)
        }
        if (isTRUE(quadrants)) {
            quadrants <- c(median(df$.x), median(df$.y))
        }
        .plot_subject_scatter(df, labels, colour.by, quadrants, ...)
    }
)

# One row per subject: .subject, .x, .y and colour.by from colData.
.subject_stats <- function(x, x.stat, y.stat, group, colour.by) {
    for (arg in c("x.stat", "y.stat")) {
        col <- get(arg)
        if (is.null(col)) {
            stop("'", arg, "' must be one column name of colData(x).",
                call. = FALSE)
        }
        .check_column(x, col, arg = arg)
        if (!is.numeric(colData(x)[[col]])) {
            stop("'", arg, "' must be a numeric column.", call. = FALSE)
        }
    }
    cols <- unique(c(group, x.stat, y.stat, colour.by))
    df <- unique(as.data.frame(colData(x)[, cols, drop = FALSE],
        optional = TRUE))
    if (anyDuplicated(df[[group]])) {
        stop("'x.stat', 'y.stat' and 'colour.by' must be constant ",
            "within each subject of 'group'.", call. = FALSE)
    }
    df$.subject <- df[[group]]
    df$.x <- df[[x.stat]]
    df$.y <- df[[y.stat]]
    attr(df, "labels") <- c(x.stat, y.stat)
    df
}

# Subject loadings of two components of a decomposition, matched to the
# subjects of group, with colour.by.
.subject_loading_table <- function(x, source, ncomponents, group,
        colour.by) {
    if (!.is_one_name(source) || !source %in% names(metadata(x))) {
        stop("'source' must name a decomposition in metadata(x).",
            call. = FALSE)
    }
    obj <- metadata(x)[[source]]
    mat <- if (is.list(obj) && !is.data.frame(obj) && "A_hat" %in%
            names(obj)) {
        as.matrix(obj$A_hat)
    } else if (is.list(obj) && !is.data.frame(obj) && "Fac" %in%
            names(obj)) {
        as.matrix(obj$Fac[[1L]])
    } else if (is.matrix(obj) || is.data.frame(obj)) {
        as.matrix(obj)
    } else {
        stop("'source' is not a TEMPTED result, a PARAFAC model or a ",
            "loading matrix.", call. = FALSE)
    }
    if (is.null(rownames(mat))) {
        stop("The subject loadings of 'source' need subjects as row ",
            "names.", call. = FALSE)
    }
    j <- .components(ncomponents, mat, n = 2L)
    cols <- unique(c(group, colour.by))
    cd <- unique(as.data.frame(colData(x)[, cols, drop = FALSE],
        optional = TRUE))
    if (anyDuplicated(cd[[group]])) {
        stop("'colour.by' must be constant within each subject of ",
            "'group'.", call. = FALSE)
    }
    i <- match(as.character(cd[[group]]), rownames(mat))
    if (all(is.na(i))) {
        stop("No subject of 'group' is in the loadings of 'source'.",
            call. = FALSE)
    }
    cd <- cd[!is.na(i), , drop = FALSE]
    i <- i[!is.na(i)]
    cd$.subject <- cd[[group]]
    cd$.x <- mat[i, j[1L]]
    cd$.y <- mat[i, j[2L]]
    attr(cd, "labels") <- .component_labels(obj, mat)[j]
    cd
}

#' @importFrom ggplot2 ggplot aes .data geom_point geom_hline geom_vline
#'   labs
.plot_subject_scatter <- function(df, labels, colour.by, quadrants, ...) {
    p <- ggplot(df, aes(x = .data$.x, y = .data$.y))
    if (is.numeric(quadrants)) {
        p <- p + geom_vline(xintercept = quadrants[[1L]], linetype = "dashed",
                colour = "grey40") +
            geom_hline(yintercept = quadrants[[2L]], linetype = "dashed",
                colour = "grey40")
    }
    p + geom_point(aes(colour = if (!is.null(colour.by))
            .data[[colour.by]]), ...) +
        labs(x = labels[[1L]], y = labels[[2L]], colour = colour.by)
}
