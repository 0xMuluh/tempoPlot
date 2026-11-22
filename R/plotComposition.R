#' Plot composition over time
#'
#' \code{plotComposition()} stacks the relative abundances of the most
#' abundant features over time, for each subject or averaged over groups of
#' subjects.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. Each subject gets its own
#'   panel. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay of
#'   abundances. To plot taxonomic ranks, agglomerate \code{x} first with
#'   \code{\link[mia:agglomerate-methods]{mia::agglomerateByRank()}}.
#'   (Default: \code{"counts"})
#'
#' @param features \code{Character vector}. Features to plot. If
#'   \code{NULL}, plots the \code{top} features with the largest mean
#'   plotted value: proportion, or assay value with
#'   \code{as.relative = FALSE}. The rest, and a feature already named
#'   "Other", are summed into "Other". (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of features to plot when
#'   \code{features} is \code{NULL}. (Default: \code{5})
#'
#' @param as.relative \code{Logical scalar}. Turns the assay into
#'   proportions of each sample (\code{TRUE}) or stacks its values as they
#'   are (\code{FALSE}), for absolute abundances such as relative
#'   abundances times a qPCR or flow cytometry load. (Default: \code{TRUE})
#'
#' @param layout \code{Character scalar}. Draws stacked bars (\code{"bar"}),
#'   stacked areas (\code{"area"}) or areas stacked around zero
#'   (\code{"stream"}), which shows changes of the total with
#'   \code{as.relative = FALSE}. (Default: \code{"bar"})
#'
#' @param level \code{Character scalar}. Plots each subject
#'   (\code{"subject"}) or the mean composition of each \code{group.by}
#'   group at each time (\code{"group"}). (Default: \code{"subject"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that groups subjects when \code{level = "group"}. If
#'   \code{NULL}, all subjects form one group. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_bar]{geom_col}},
#'   \code{\link[ggplot2:geom_ribbon]{geom_area}} or
#'   \code{\link[ggplot2:geom_ribbon]{geom_ribbon}}.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotComposition
#'
#' @examples
#' # A defined community of 16 species in three replicate bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- minimalgut
#' plotComposition(tse, time.col = "Time_hr_num",
#'     group = "StudyIdentifier", top = 8, layout = "area")
#'
#' # Biomass of each species, relative abundance times optical density,
#' # as a stream
#' tse <- tse[, !is.na(tse$OpticalDensity)]
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' SummarizedExperiment::assay(tse, "biomass") <-
#'     t(t(SummarizedExperiment::assay(tse, "relabundance")) *
#'     tse$OpticalDensity)
#' plotComposition(tse, time.col = "Time_hr_num",
#'     group = "StudyIdentifier", assay.type = "biomass", top = 8,
#'     as.relative = FALSE, layout = "stream")
#'
#' # Phyla of five subjects over their visits
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' tse <- mia::agglomerateByRank(tse, rank = "Phylum")
#' plotComposition(tse, time.col = "visit", group = "subject")
#'
#' # Mean phylum composition by sex at each visit
#' plotComposition(tse, time.col = "visit", group = "subject",
#'     level = "group", group.by = "sex")
NULL

#' @rdname plotComposition
#' @importFrom SummarizedExperiment colData
#' @importFrom stats aggregate
#' @importFrom methods setMethod
#' @export
setMethod("plotComposition", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            features = NULL, top = 5, as.relative = TRUE,
            layout = c("bar", "area", "stream"),
            level = c("subject", "group"), group.by = NULL, ...) {
        layout <- match.arg(layout)
        level <- match.arg(level)
        .check_longitudinal(x, time.col, group)
        .check_column(x, group.by)
        if (!(isTRUE(as.relative) || isFALSE(as.relative))) {
            stop("'as.relative' must be TRUE or FALSE.", call. = FALSE)
        }
        mat <- .get_assay(x, assay.type)
        if (as.relative) {
            mat <- prop.table(mat, 2L)
        }
        mat <- .sum_other(mat, features, top)

        cols <- unique(c(time.col, group, group.by))
        cd <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df <- cd[rep(seq_len(ncol(x)), each = nrow(mat)), , drop = FALSE]
        df$.feature <- factor(rownames(mat), levels = rownames(mat))
        df$.value <- as.vector(mat)
        panel <- group
        if (level == "group") {
            df <- aggregate(df[".value"], df[c(time.col, group.by,
                ".feature")], mean)
            panel <- group.by
        }
        ylab <- if (as.relative) "Relative abundance" else assay.type
        .plot_composition(df, time.col, panel, layout, ylab, ...)
    }
)

# Rows of the selected features, in selection order, and a last row
# "Other" with the sum of the rest. A row already named "Other" is never
# selected and is summed with the rest.
.sum_other <- function(mat, features, top) {
    keep <- .select_features(mat[rownames(mat) != "Other", , drop = FALSE],
        features, top)
    lab <- rownames(mat)
    lab[!lab %in% keep] <- "Other"
    rowsum(mat, factor(lab, levels = c(keep, "Other")))
}

#' @importFrom ggplot2 ggplot aes .data geom_col geom_area geom_ribbon
#'   facet_wrap vars labs scale_x_discrete scale_fill_manual theme
#'   element_blank
#' @importFrom stats ave
.plot_composition <- function(df, time.col, panel, layout, ylab, ...) {
    if (layout == "stream") {
        # Stack the features at each time, the first on top as in bars and
        # areas, and centre the stack on zero. Only heights are meaningful,
        # so the y axis is left unlabelled.
        df <- df[order(df[[time.col]], -as.integer(df$.feature)), ,
            drop = FALSE]
        keys <- if (is.null(panel)) list(df[[time.col]]) else
            list(df[[time.col]], df[[panel]])
        top <- ave(df$.value, keys, FUN = cumsum)
        total <- ave(df$.value, keys, FUN = sum)
        df$.ymax <- top - total / 2
        df$.ymin <- df$.ymax - df$.value
        p <- ggplot(df, aes(x = .data[[time.col]], ymin = .data$.ymin,
                ymax = .data$.ymax, fill = .data$.feature,
                group = .data$.feature)) +
            geom_ribbon(...)
    } else {
        geom <- switch(layout, bar = geom_col, area = geom_area)
        p <- ggplot(df, aes(x = .data[[time.col]], y = .data$.value,
                fill = .data$.feature, group = .data$.feature)) +
            geom(...)
    }
    # The features in the categorical palette, "Other" in neutral grey.
    lv <- setdiff(levels(df$.feature), "Other")
    values <- c(.palette_categorical(length(lv)), .na_colour)
    names(values) <- c(lv, "Other")
    p <- p + scale_fill_manual(values = values) +
        labs(x = time.col, y = ylab, fill = "feature") +
        .theme_tempoPlot()
    if (layout == "stream") {
        p <- p + theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    }
    if (!is.null(panel)) {
        p <- p + facet_wrap(vars(.data[[panel]]))
    }
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p
}
