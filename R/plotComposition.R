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
#' @param layout \code{Character scalar}. Draws stacked bars (\code{"bar"})
#'   or stacked areas (\code{"area"}). (Default: \code{"bar"})
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
#'   \code{\link[ggplot2:geom_bar]{geom_col}} or
#'   \code{\link[ggplot2:geom_ribbon]{geom_area}}.
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
NULL

#' @rdname plotComposition
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom stats aggregate
#' @importFrom methods setMethod
#' @export
setMethod("plotComposition", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = NULL, rank = NULL,
            top = 5, layout = c("bar", "area"),
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
        if (level == "group" && is.null(group.by)) {
            stop("'level' = \"group\" needs 'group.by'.", call. = FALSE)
        }

        mat <- as.matrix(assay(x, assay.type %||% 1L))
        if (!is.null(rank)) {
            mat <- rowsum(mat, as.character(rowData(x)[[rank]]))
        }
        prop <- prop.table(mat, 2L)
        keep <- order(rowMeans(prop), decreasing = TRUE)[seq_len(
            min(top, nrow(prop)))]
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
            rank %||% "feature", ...)
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_col geom_area facet_wrap vars
#'   labs scale_x_discrete
.plot_composition <- function(df, time.col, panel, layout, fill, ...) {
    geom <- switch(layout, bar = geom_col, area = geom_area)
    p <- ggplot(df, aes(x = .data[[time.col]], y = .data$.value,
            fill = .data$.feature, group = .data$.feature)) +
        geom(...) +
        labs(x = time.col, y = "Relative abundance", fill = fill)
    if (!is.null(panel)) {
        p <- p + facet_wrap(vars(.data[[panel]]))
    }
    if (is.factor(df[[time.col]])) {
        p <- p + scale_x_discrete(drop = FALSE)
    }
    p
}
