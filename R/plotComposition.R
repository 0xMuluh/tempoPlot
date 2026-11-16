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
#'   abundances. Values are turned into proportions of each sample.
#'   (Default: \code{NULL}, the first assay)
#'
#' @param rank \code{Character scalar}. Specifies the column of
#'   \code{rowData(x)}, such as a taxonomic rank, whose levels are summed
#'   before plotting. (Default: \code{NULL})
#'
#' @param top \code{Integer scalar}. Number of most abundant features, by
#'   mean proportion, to show; the rest are summed into "Other".
#'   (Default: \code{5})
#'
#' @param y \code{Character scalar}. Specifies a categorical column from
#'   \code{colData(x)}, such as a state, to plot instead of the assay: the
#'   share of samples in each state at each time, one panel per
#'   \code{group.by} group. (Default: \code{NULL})
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
#'   \code{colData(x)} that groups subjects when \code{level = "group"}.
#'   (Default: \code{NULL})
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
#' plotComposition(minimalgut, time.col = "Time_hr_num",
#'     group = "StudyIdentifier", top = 8, layout = "area")
#'
#' # Phyla of five subjects over their visits
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' plotComposition(tse, time.col = "visit", group = "subject",
#'     rank = "Phylum")
#'
#' # Mean phylum composition by sex at each visit
#' plotComposition(tse, time.col = "visit", group = "subject",
#'     rank = "Phylum", level = "group", group.by = "sex")
#'
#' # Biomass of each species, relative abundance times optical density,
#' # as a stream
#' mg <- minimalgut[, !is.na(minimalgut$OpticalDensity)]
#' mg <- mia::transformAssay(mg, method = "relabundance")
#' SummarizedExperiment::assay(mg, "biomass") <-
#'     t(t(SummarizedExperiment::assay(mg, "relabundance")) *
#'     mg$OpticalDensity)
#' plotComposition(mg, time.col = "Time_hr_num",
#'     group = "StudyIdentifier", assay.type = "biomass", top = 8,
#'     as.relative = FALSE, layout = "stream")
#'
#' # Share of samples by dominant genus at each visit
#' tse <- mia::addDominant(tse, name = "dominant", n = 4)
#' plotComposition(tse, time.col = "visit", group = "subject",
#'     y = "dominant")
NULL

#' @rdname plotComposition
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom stats aggregate
#' @importFrom methods setMethod
#' @export
setMethod("plotComposition", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = NULL, rank = NULL,
            top = 5, y = NULL, as.relative = TRUE,
            layout = c("bar", "area", "stream"),
            level = c("subject", "group"), group.by = NULL, ...) {
        layout <- match.arg(layout)
        level <- match.arg(level)
        .check_longitudinal(x, time.col, group)
        # Errors on character and other unsupported time types.
        .time_type(x, time.col)
        .check_column(x, group.by)
        .check_column(x, rank, where = "rowData")
        if (!(is.numeric(top) && length(top) == 1L && top >= 1)) {
            stop("'top' must be a single number of at least 1.",
                call. = FALSE)
        }
        if (!(isTRUE(as.relative) || isFALSE(as.relative))) {
            stop("'as.relative' must be TRUE or FALSE.", call. = FALSE)
        }
        if (level == "group" && is.null(group.by)) {
            stop("'level' = \"group\" needs 'group.by'.", call. = FALSE)
        }
        if (!is.null(y)) {
            df <- .state_shares(x, y, time.col, group.by)
            return(.plot_composition(df, time.col, group.by, layout, y,
                ylab = "Share of samples", ...))
        }

        mat <- as.matrix(assay(x, assay.type %||% 1L))
        if (!is.null(rank)) {
            mat <- rowsum(mat, as.character(rowData(x)[[rank]]))
        }
        prop <- prop.table(mat, 2L)
        keep <- order(rowMeans(prop), decreasing = TRUE)[seq_len(
            min(top, nrow(prop)))]
        if (!as.relative) {
            prop <- mat
        }
        if (length(keep) < nrow(prop)) {
            prop <- rbind(prop[keep, , drop = FALSE],
                Other = colSums(prop[-keep, , drop = FALSE]))
        }

        cols <- unique(c(time.col, group, group.by))
        cd <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        df <- cd[rep(seq_len(ncol(x)), each = nrow(prop)), , drop = FALSE]
        df$.feature <- factor(rownames(prop), levels = rownames(prop))
        df$.value <- as.vector(prop)
        panel <- group
        if (level == "group") {
            df <- aggregate(df[".value"], df[c(time.col, group.by,
                ".feature")], mean)
            panel <- group.by
        }

        .plot_composition(df, time.col, panel, layout,
            rank %||% "feature", ylab = if (as.relative)
                "Relative abundance" else assay.type %||% "Abundance",
            ...)
    }
)

# Share of samples in each state of the colData column y at each time,
# within each group.by group.
#' @importFrom SummarizedExperiment colData
.state_shares <- function(x, y, time.col, group.by) {
    .check_column(x, y)
    if (is.null(y)) {
        stop("'y' must be one column name of colData(x).", call. = FALSE)
    }
    state <- colData(x)[[y]]
    if (!(is.factor(state) || is.character(state) || is.logical(state))) {
        stop("'y' must be a factor, character or logical column.",
            call. = FALSE)
    }
    cols <- unique(c(time.col, group.by))
    cd <- as.data.frame(colData(x)[, cols, drop = FALSE], optional = TRUE)
    cd$.feature <- factor(state)
    cd <- cd[!is.na(cd$.feature), , drop = FALSE]
    tab <- as.data.frame(table(cd[c(time.col, group.by, ".feature")]),
        responseName = ".n")
    tab$.value <- tab$.n / ave(tab$.n, tab[c(time.col, group.by)],
        FUN = sum)
    tab <- tab[!is.nan(tab$.value), , drop = FALSE]
    # table() turns times into factors; restore numeric and Date times.
    time <- colData(x)[[time.col]]
    if (!is.factor(time)) {
        lv <- as.character(tab[[time.col]])
        tab[[time.col]] <- sort(unique(time))[match(lv,
            as.character(sort(unique(time))))]
    }
    tab
}

#' @importFrom ggplot2 ggplot aes .data geom_col geom_area geom_ribbon
#'   facet_wrap vars labs scale_x_discrete
#' @importFrom stats ave
.plot_composition <- function(df, time.col, panel, layout, fill,
        ylab = "Relative abundance", ...) {
    if (layout == "stream") {
        # Stack the features at each time, then centre the stack on zero.
        df <- df[order(df[[time.col]], df$.feature), , drop = FALSE]
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
    p <- p + labs(x = time.col, y = ylab, fill = fill)
    if (!is.null(panel)) {
        p <- p + facet_wrap(vars(.data[[panel]]))
    }
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p
}
