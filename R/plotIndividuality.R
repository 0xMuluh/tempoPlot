#' Plot within- against between-subject dissimilarities
#'
#' \code{plotIndividuality()} compares the dissimilarities between samples
#' of the same subject with those between samples of different subjects.
#' The further apart the two distributions, the more individual the
#' communities.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param dissimilarity \code{Character scalar}. Specifies the
#'   sample by sample dissimilarity in \code{metadata(x)}, as stored by
#'   \code{\link[mia:getDissimilarity]{mia::addDissimilarity()}}.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param layout \code{Character scalar}. Draws density curves
#'   (\code{"density"}), box plots (\code{"box"}), violins
#'   (\code{"violin"}) or the median and interquartile range in each time
#'   window (\code{"time"}). (Default: \code{"density"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the between-subject pairs into pairs of
#'   the same and of different groups, such as the study centre.
#'   (Default: \code{NULL})
#'
#' @param time.col \code{Character scalar}. For \code{layout = "time"},
#'   specifies the column from \code{colData(x)} holding the time window
#'   of each sample: a factor, or numeric, \code{Date} or \code{POSIXct}
#'   values of which each is one window. (Default: \code{NULL})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that splits the pairs into panels. Only pairs of
#'   samples from the same group are kept. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_density]{geom_density}},
#'   \code{\link[ggplot2:geom_boxplot]{geom_boxplot}},
#'   \code{\link[ggplot2:geom_violin]{geom_violin}} or, for
#'   \code{layout = "time"}, \code{\link[ggplot2:geom_line]{geom_line}}.
#'
#' @details
#' Each pair of samples is counted once. A subject with \eqn{k} samples
#' gives \eqn{k(k-1)/2} within-subject pairs, so subjects sampled often
#' weigh most in the within distribution.
#'
#' With \code{layout = "time"}, only pairs of samples from the same window
#' are kept, and each line is the median dissimilarity of a kind of pair
#' across the windows, as in Vatanen et al. (2018,
#' doi:10.1038/s41586-018-0620-2, Fig. 2b). Make windows from times with
#' \code{\link[base:cut]{cut()}}.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotIndividuality
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 2))]
#' tse <- mia::addDissimilarity(tse, method = "bray", assay.type = "counts")
#' plotIndividuality(tse, dissimilarity = "bray", group = "subject")
#'
#' # By sex
#' plotIndividuality(tse, dissimilarity = "bray", group = "subject",
#'     layout = "box", facet.by = "sex")
#'
#' # Across age windows, with pairs of subjects of the same and of
#' # different sex
#' tse$age.window <- cut(tse$age, c(18, 30, 40, 50, 60, 80),
#'     include.lowest = TRUE)
#' plotIndividuality(tse, dissimilarity = "bray", group = "subject",
#'     layout = "time", time.col = "age.window", group.by = "sex")
NULL

#' @rdname plotIndividuality
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotIndividuality", signature = c(x = "SummarizedExperiment"),
    function(x, dissimilarity, group,
            layout = c("density", "box", "violin", "time"), group.by = NULL,
            time.col = NULL, facet.by = NULL, ...) {
        layout <- match.arg(layout)
        .require_column(x, group)
        .check_column(x, group.by)
        .check_column(x, facet.by)
        if (layout == "time") {
            .require_column(x, time.col)
            .time_type(colData(x)[[time.col]])
        } else {
            time.col <- NULL
        }
        d <- .get_dissimilarity(x, dissimilarity)

        pairs <- which(lower.tri(d), arr.ind = TRUE)
        df <- data.frame(dissimilarity = d[pairs],
            pair = .pair_type(colData(x), pairs, group, group.by))
        # Keep pairs whose two samples share the window and the panel.
        keep <- !is.na(df$dissimilarity) & !is.na(df$pair)
        for (col in c(time.col, facet.by)) {
            v <- colData(x)[[col]]
            same <- v[pairs[, 1L]] == v[pairs[, 2L]]
            keep <- keep & !is.na(same) & same
            df[[col]] <- v[pairs[, 1L]]
        }
        df <- df[keep, , drop = FALSE]
        if (!any(df$pair == "within")) {
            stop("No subject has more than one sample.", call. = FALSE)
        }
        .plot_individuality(df, dissimilarity, layout, time.col, facet.by,
            ...)
    }
)

# Kind of each pair of samples (rows of pairs): within a subject, or
# between subjects, split by whether they share group.by. NA when
# group.by is missing for a between-subject pair.
.pair_type <- function(cd, pairs, group, group.by) {
    same <- function(col) {
        v <- cd[[col]]
        v[pairs[, 1L]] == v[pairs[, 2L]]
    }
    within <- same(group)
    if (is.null(group.by)) {
        return(factor(ifelse(within, "within", "between"),
            levels = c("within", "between")))
    }
    levels <- c("within", paste(c("between, same", "between, different"),
        group.by))
    factor(ifelse(within, levels[[1L]],
        ifelse(same(group.by), levels[[2L]], levels[[3L]])), levels = levels)
}

#' @importFrom ggplot2 ggplot aes .data geom_density geom_boxplot
#'   geom_violin labs facet_wrap vars
.plot_individuality <- function(df, dissimilarity, layout, time.col,
        facet.by, ...) {
    if (layout == "density") {
        p <- ggplot(df, aes(x = .data$dissimilarity,
                colour = .data$pair)) +
            geom_density(...) +
            .scale_categorical("colour") +
            labs(x = dissimilarity, colour = "pair")
    } else if (layout == "time") {
        p <- ggplot() +
            .fit_layers(.pair_quartiles(df, time.col, facet.by), time.col,
                "pair", ...) +
            .scale_categorical(c("colour", "fill")) +
            labs(x = time.col, y = dissimilarity, colour = "pair",
                fill = "pair")
    } else {
        geom <- if (layout == "box") geom_boxplot else geom_violin
        p <- ggplot(df, aes(x = .data$pair, y = .data$dissimilarity)) +
            geom(...) +
            labs(x = "pair", y = dissimilarity)
    }
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data[[facet.by]]))
    }
    p + .theme_tempoPlot()
}

# Median (estimate) and quartiles (lower, upper) of the dissimilarities
# of each kind of pair in each window and panel.
#' @importFrom stats aggregate quantile
.pair_quartiles <- function(df, time.col, facet.by) {
    q <- aggregate(df["dissimilarity"], df[c(time.col, "pair", facet.by)],
        quantile, probs = c(0.25, 0.5, 0.75), names = FALSE)
    q[c("lower", "estimate", "upper")] <- as.data.frame(q$dissimilarity)
    q$dissimilarity <- NULL
    q
}
