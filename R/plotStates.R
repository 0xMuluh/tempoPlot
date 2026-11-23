#' Plot the share of samples of each subject in each state
#'
#' \code{plotStates()} draws, for each subject, the share of its samples in
#' each state, such as a community state type, as a subject by state
#' heatmap or as stacked bars.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param state.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains the state of each sample, such as a
#'   community type or the dominant taxon.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)}, constant within subjects, that splits the subjects
#'   into panels. (Default: \code{NULL})
#'
#' @param layout \code{Character scalar}. Draws a subject by state heatmap
#'   (\code{"heatmap"}) or one stacked bar per subject (\code{"bar"}).
#'   (Default: \code{"heatmap"})
#'
#' @param order.by \code{Character scalar}. Orders subjects by their most
#'   common state, then by its share (\code{"dominant"}), or by
#'   hierarchical clustering of their shares (\code{"cluster"}).
#'   (Default: \code{"dominant"})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_tile]{geom_tile}} or
#'   \code{\link[ggplot2:geom_col]{geom_col}}.
#'
#' @details
#' Shares are of the samples with a state: samples without one are left
#' out. Each sample counts once, whatever the time to the next sample, so
#' subjects should be sampled at regular times. The heatmap follows the
#' community state type frequencies of Gajer et al. (2012), clustered into
#' community classes, and the bars those of DiGiulio et al. (2015), split
#' by pregnancy outcome.
#'
#' @references
#' Gajer P, et al. (2012). Temporal dynamics of the human vaginal
#' microbiota. \emph{Science Translational Medicine} 4(132):132ra52.
#' doi:10.1126/scitranslmed.3003605
#'
#' DiGiulio DB, et al. (2015). Temporal and spatial variation of the human
#' microbiota during pregnancy. \emph{PNAS} 112(35):11060-11065.
#' doi:10.1073/pnas.1502875112
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotStates
#'
#' @examples
#' # Share of daily samples of each person in each enterotype. Samples
#' # without a cell count have the enterotype "NA".
#' data("temporalMicrobiome20", package = "miaTime")
#' tse <- temporalMicrobiome20
#' tse$enterotype[tse$enterotype == "NA"] <- NA
#' plotStates(tse, state.col = "enterotype", group = "person_id")
#'
#' # Dominant genus of the subjects with at least five samples, by sex
#' data("hitchip1006", package = "miaTime")
#' tse <- mia::addDominant(hitchip1006, name = "dominant", n = 4)
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 5))]
#' plotStates(tse, state.col = "dominant", group = "subject",
#'     group.by = "sex", layout = "bar", order.by = "cluster")
NULL

#' @rdname plotStates
#' @importFrom SummarizedExperiment colData
#' @importFrom methods setMethod
#' @export
setMethod("plotStates", signature = c(x = "SummarizedExperiment"),
    function(x, state.col, group, group.by = NULL,
            layout = c("heatmap", "bar"),
            order.by = c("dominant", "cluster"), ...) {
        layout <- match.arg(layout)
        order.by <- match.arg(order.by)
        .require_column(x, state.col)
        .require_column(x, group)
        .check_column(x, group.by)

        state <- colData(x)[[state.col]]
        keep <- !is.na(state)
        if (!any(keep)) {
            stop("No sample has a state.", call. = FALSE)
        }
        subject <- as.character(colData(x)[[group]])[keep]
        share <- prop.table(table(subject, factor(as.character(state[keep]),
            levels = .levels_of(state))), 1L)
        df <- as.data.frame(share)
        names(df) <- c(".subject", ".state", ".share")
        # The first subject at the top.
        df$.subject <- factor(df$.subject,
            levels = rev(.order_subjects(share, order.by)))
        if (!is.null(group.by)) {
            by <- as.character(colData(x)[[group.by]])[keep]
            pairs <- unique(data.frame(subject, by))
            if (anyDuplicated(pairs$subject)) {
                stop("'group.by' must be constant within subjects.",
                    call. = FALSE)
            }
            df$.by <- pairs$by[match(df$.subject, pairs$subject)]
        }
        .plot_states(df, state.col, group, group.by, layout, ...)
    }
)

# Rows of the subject by state share matrix, in the order of order.by.
#' @importFrom stats dist hclust
.order_subjects <- function(share, order.by) {
    if (order.by == "cluster") {
        if (nrow(share) < 2L) {
            return(rownames(share))
        }
        return(rownames(share)[hclust(dist(share))$order])
    }
    top <- max.col(share, ties.method = "first")
    rownames(share)[order(top, -share[cbind(seq_along(top), top)])]
}

#' @importFrom ggplot2 ggplot aes .data geom_tile geom_col labs facet_grid
#'   vars theme element_text scale_x_continuous
.plot_states <- function(df, state.col, group, group.by, layout, ...) {
    if (layout == "heatmap") {
        p <- ggplot(df, aes(x = .data$.state, y = .data$.subject)) +
            geom_tile(aes(fill = .data$.share), ...) +
            .scale_sequential("fill") +
            labs(x = state.col, fill = "share") + .theme_tempoPlot() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
        p <- ggplot(df, aes(x = .data$.share, y = .data$.subject)) +
            geom_col(aes(fill = .data$.state), ...) +
            .scale_categorical("fill") +
            scale_x_continuous(expand = c(0, 0)) +
            labs(x = "share of samples", fill = state.col) +
            .theme_tempoPlot()
    }
    p <- p + labs(y = group)
    if (!is.null(group.by)) {
        p <- p + facet_grid(rows = vars(.data$.by), scales = "free_y",
            space = "free_y")
    }
    p
}
