#' Plot lagged correlations of features with a target
#'
#' \code{plotLagHeatmap()} correlates each feature with a target series
#' shifted by a number of samples within each subject, and draws the
#' correlations as features by lag tiles.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   or \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object.
#'   For a \code{MultiAssayExperiment}, \code{time.col} and \code{group}
#'   are columns of its \code{colData}.
#'
#' @param target \code{Character scalar}. Specifies the target series: a
#'   column of \code{colData(x)}, a row of \code{assay(x, assay.type)}, or
#'   \code{"<dimred>:<component>"}.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of the
#'   features. (Default: \code{"counts"})
#'
#' @param features \code{Character vector}. Features to correlate. If
#'   \code{NULL}, uses the \code{top} features with the largest mean.
#'   (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features when
#'   \code{features} is \code{NULL}. (Default: \code{20})
#'
#' @param experiment,target.experiment \code{Character scalar}. For a
#'   \code{MultiAssayExperiment}, the experiments holding the features and
#'   the \code{target}; samples are matched by their primary names. The
#'   target is looked up in \code{target.experiment} as for a
#'   \code{SummarizedExperiment}, in its first assay for a feature.
#'
#' @param max.lag \code{Integer scalar}. Largest lag in samples, in both
#'   directions. (Default: \code{3})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_tile]{geom_tile}}.
#'
#' @details
#' At lag \eqn{k}, each feature at a sample is paired with the target
#' \eqn{k} samples later in the same subject, and the Spearman correlation
#' is computed over all pairs. Positive lags show features that change
#' before the target, negative lags features that follow it. Lags count
#' samples, so they are comparable across subjects only for even sampling.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotLagHeatmap
#'
#' @examples
#' # Species that change before or after acetate in three bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- mia::transformAssay(minimalgut, method = "clr", pseudocount = 1)
#' plotLagHeatmap(tse, target = "Acetate", time.col = "Time.hr",
#'     group = "StudyIdentifier", assay.type = "clr")
#'
#' # The same with acetate from a separate experiment
#' cd <- SummarizedExperiment::colData(minimalgut)
#' scfa <- SummarizedExperiment::SummarizedExperiment(
#'     t(as.matrix(cd[, c("Acetate", "Butyrate", "Propionate")])))
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'     experiments = list(species = tse, SCFA = scfa), colData = cd)
#' plotLagHeatmap(mae, target = "Acetate", time.col = "Time.hr",
#'     group = "StudyIdentifier", experiment = "species",
#'     target.experiment = "SCFA", assay.type = "clr")
NULL

#' @rdname plotLagHeatmap
#' @importFrom SummarizedExperiment colData assay
#' @importFrom methods setMethod
#' @export
setMethod("plotLagHeatmap", signature = c(x = "SummarizedExperiment"),
    function(x, target, time.col, group = NULL, assay.type = "counts",
            features = NULL, top = 20, max.lag = 3, ...) {
        df <- .feature_table(x, time.col, group, assay.type, features, top)
        y <- unname(.resolve_y(x, target, NULL, "target"))
        mat <- as.matrix(assay(x, assay.type))[unique(df$.feature), ,
            drop = FALSE]
        time <- .time_numeric(colData(x)[[time.col]])
        subject <- if (is.null(group)) rep("", ncol(x)) else
            as.character(colData(x)[[group]])
        res <- .lag_correlations(mat, y, time, subject, max.lag)
        .plot_lag_heatmap(res, target, ...)
    }
)

#' @rdname plotLagHeatmap
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom methods setMethod
#' @export
setMethod("plotLagHeatmap", signature = c(x = "MultiAssayExperiment"),
    function(x, target, time.col, group = NULL, experiment,
            target.experiment, assay.type = "counts", features = NULL,
            top = 20, max.lag = 3, ...) {
        exps <- names(experiments(x))
        for (arg in c("experiment", "target.experiment")) {
            if (!.is_one_name(get(arg)) || !get(arg) %in% exps) {
                stop("'", arg, "' must be one experiment of 'x'.",
                    call. = FALSE)
            }
        }
        cd <- colData(x)
        for (arg in c("time.col", "group")) {
            v <- get(arg)
            if (arg == "group" && is.null(v)) next
            if (!.is_one_name(v) || !v %in% names(cd)) {
                stop("'", arg, "' must be one column name of colData(x).",
                    call. = FALSE)
            }
        }
        e1 <- x[[experiment]]
        e2 <- x[[target.experiment]]
        if (!.is_one_name(assay.type) || !assay.type %in% assayNames(e1)) {
            stop("'assay.type' must be one of the assays of 'experiment'.",
                call. = FALSE)
        }
        mat <- as.matrix(assay(e1, assay.type))
        mat <- mat[.select_features(mat, features, top), , drop = FALSE]
        y <- .resolve_y(e2, target, NULL, "target")

        # Match the samples of both experiments by their primary names.
        map <- as.data.frame(sampleMap(x))
        primary <- function(e, cols) {
            m <- map[map$assay == e, , drop = FALSE]
            m$primary[match(cols, m$colname)]
        }
        p1 <- primary(experiment, colnames(e1))
        p2 <- primary(target.experiment, colnames(e2))
        common <- intersect(p1, p2)
        if (anyDuplicated(p1[p1 %in% common]) ||
                anyDuplicated(p2[p2 %in% common])) {
            stop("Each primary sample must have one sample in each ",
                "experiment.", call. = FALSE)
        }
        if (length(common) < 3L) {
            stop("The experiments share fewer than three samples.",
                call. = FALSE)
        }
        mat <- mat[, match(common, p1), drop = FALSE]
        y <- unname(y[match(common, p2)])
        time <- cd[common, time.col]
        if (!(is.numeric(time) || inherits(time, "Date") ||
                is.factor(time))) {
            stop("'time.col' must be numeric, Date or factor.",
                call. = FALSE)
        }
        subject <- if (is.null(group)) rep("", length(common)) else
            as.character(cd[common, group])
        if (anyDuplicated(data.frame(subject, time))) {
            stop("'time.col' must not repeat within a subject.",
                call. = FALSE)
        }
        res <- .lag_correlations(mat, y, .time_numeric(time), subject,
            max.lag)
        .plot_lag_heatmap(res, target, ...)
    }
)

# Spearman correlation of each feature (rows of mat) with y at lags of
# -max.lag to max.lag samples within each subject. Positive lags pair a
# feature with y later in the subject.
#' @importFrom stats cor
.lag_correlations <- function(mat, y, time, subject, max.lag) {
    if (!(is.numeric(max.lag) && length(max.lag) == 1L && max.lag >= 0)) {
        stop("'max.lag' must be a single non-negative number.",
            call. = FALSE)
    }
    max.lag <- as.integer(max.lag)
    ord <- order(subject, time)
    n <- length(ord)
    res <- do.call(rbind, lapply(-max.lag:max.lag, function(k) {
        a <- seq_len(n)
        b <- a + k
        ok <- b >= 1L & b <= n
        a <- a[ok]
        b <- b[ok]
        same <- subject[ord[a]] == subject[ord[b]]
        a <- ord[a[same]]
        b <- ord[b[same]]
        if (length(a) < 3L) return(NULL)
        r <- apply(mat[, a, drop = FALSE], 1L, function(f)
            suppressWarnings(cor(f, y[b], method = "spearman",
                use = "pairwise.complete.obs")))
        data.frame(.feature = rownames(mat), .lag = k, .r = r)
    }))
    if (is.null(res)) {
        stop("Too few pairs of samples for the lags.", call. = FALSE)
    }
    res$.feature <- factor(res$.feature, levels = rev(rownames(mat)))
    res
}

#' @importFrom ggplot2 ggplot aes .data geom_tile labs scale_fill_gradient2
.plot_lag_heatmap <- function(df, target, ...) {
    ggplot(df, aes(x = .data$.lag, y = .data$.feature, fill = .data$.r)) +
        geom_tile(...) +
        scale_fill_gradient2(low = "#b2182b", mid = "grey95",
            high = "#2166ac", limits = c(-1, 1), na.value = "white") +
        labs(x = paste0("lag in samples (feature before ", target, ")"),
            y = NULL, fill = "Spearman")
}
