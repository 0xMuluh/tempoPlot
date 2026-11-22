#' Plot gains and losses of features over time
#'
#' \code{plotTurnover()} counts the features gained and lost between
#' consecutive samples of each subject and draws them as bars spanning the
#' interval between the samples.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits. Time differences, as in
#'   \code{max.gap}, are in days for \code{Date} and hours for
#'   \code{POSIXct} times.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay.
#'   (Default: \code{"counts"})
#'
#' @param detection \code{Numeric scalar}. A feature is present in a
#'   sample when its value is above this. (Default: \code{0})
#'
#' @param type \code{Character scalar}. Draws the numbers of features
#'   gained (up) and lost (down) (\code{"gainloss"}), or turnover as the
#'   share of features present in either sample that appeared and
#'   disappeared (\code{"components"}). (Default: \code{"gainloss"})
#'
#' @param level \code{Character scalar}. Draws one panel per subject
#'   (\code{"subject"}) or averages the subjects of each \code{group.by}
#'   group at each time (\code{"group"}). (Default: \code{"subject"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that groups subjects when \code{level = "group"}. If
#'   \code{NULL}, all subjects form one group. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Skips the intervals between
#'   consecutive samples of a subject that are further apart than this, in
#'   \code{time.col} units. For visits, \code{NULL} skips intervals across
#'   visits without a sample; otherwise \code{NULL} keeps all of them.
#'   (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_rect]{geom_rect}}.
#'
#' @details
#' With \code{type = "components"}, the bars stack appearance and
#' disappearance; their sum is the total turnover of Hallett et al. (2016).
#'
#' @return
#' A \code{ggplot} object.
#'
#' @references
#' Hallett LM, et al. (2016). codyn: An R package of community dynamics
#' metrics. \emph{Methods in Ecology and Evolution} 7(10):1146-1151.
#'
#' @name plotTurnover
#'
#' @examples
#' # 16 species in three bioreactors
#' data("minimalgut", package = "miaTime")
#' tse <- minimalgut
#' plotTurnover(tse, time.col = "Time.hr", group = "StudyIdentifier")
#'
#' # Turnover components, averaged over the bioreactors
#' plotTurnover(tse, time.col = "Time.hr", group = "StudyIdentifier",
#'     type = "components", level = "group")
NULL

#' @rdname plotTurnover
#' @importFrom SummarizedExperiment colData
#' @importFrom stats aggregate
#' @importFrom methods setMethod
#' @export
setMethod("plotTurnover", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            detection = 0, type = c("gainloss", "components"),
            level = c("subject", "group"), group.by = NULL, max.gap = NULL,
            ...) {
        type <- match.arg(type)
        level <- match.arg(level)
        .check_longitudinal(x, time.col, group)
        .check_column(x, group.by)
        present <- .present(.get_assay(x, assay.type), detection)

        cols <- unique(c(time.col, group, group.by))
        cd <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        cd$.i <- seq_len(ncol(x))
        cd <- .order_time(cd, time.col, group)
        cd <- .break_gaps(cd, time.col, group, max.gap)
        subject <- if (is.null(group)) rep("", nrow(cd)) else
            as.character(cd[[group]])
        panel <- if (level == "subject") subject else if
            (is.null(group.by)) rep("", nrow(cd)) else
            as.character(cd[[group.by]])
        # Consecutive samples of a subject, not across a gap.
        i <- .next_in_subject(subject, segment = cd$.segment)
        if (length(i) == 0L) {
            stop("No subject has two consecutive samples.", call. = FALSE)
        }
        df <- .turnover_table(present[, cd$.i, drop = FALSE], i,
            .time_numeric(cd[[time.col]]), panel, type)
        if (level == "group") {
            df <- aggregate(cbind(.from, .value) ~ .panel + .time + .part,
                data = df, FUN = mean)
        }
        time <- cd[[time.col]]
        df <- .stack_turnover(df, type, time)
        panels <- level == "group" && !is.null(group.by) ||
            level == "subject" && !is.null(group)
        .plot_turnover(df, time.col, type, panels, levels(time), ...)
    }
)

# Features gained and lost between samples i and i + 1 (columns of present
# in time order): counts, or shares of the features present in either
# sample. One row per interval and part, with the panel of the later
# sample.
.turnover_table <- function(present, i, time, panel, type) {
    before <- present[, i, drop = FALSE]
    after <- present[, i + 1L, drop = FALSE]
    gained <- colSums(!before & after, na.rm = TRUE)
    lost <- colSums(before & !after, na.rm = TRUE)
    if (type == "components") {
        total <- colSums(before | after, na.rm = TRUE)
        gained <- ifelse(total > 0, gained / total, 0)
        lost <- ifelse(total > 0, lost / total, 0)
    }
    parts <- if (type == "gainloss") c("gained", "lost") else
        c("appearance", "disappearance")
    data.frame(.panel = panel[i + 1L], .from = time[i], .time = time[i + 1L],
        .part = rep(parts, each = length(i)),
        .value = c(gained, if (type == "gainloss") -lost else lost))
}

# Bar extents: gains and appearance from zero up; losses from zero down and
# disappearance stacked on appearance. Visit bars leave a small space
# between intervals; other times are returned in the class of time.
#' @importFrom stats ave
.stack_turnover <- function(df, type, time) {
    df <- df[order(df$.panel, df$.time, df$.part), , drop = FALSE]
    df$.ymin <- pmin(df$.value, 0)
    df$.ymax <- pmax(df$.value, 0)
    if (type == "components") {
        df$.ymin <- ave(df$.value, df$.panel, df$.time, FUN = function(v)
            c(0, cumsum(v)[-length(v)]))
        df$.ymax <- df$.ymin + df$.value
    }
    if (is.factor(time)) {
        df$.from <- df$.from + 0.05
        df$.time <- df$.time - 0.05
    }
    df$.from <- .time_restore(df$.from, time)
    df$.time <- .time_restore(df$.time, time)
    df
}

#' @importFrom utils modifyList
#' @importFrom ggplot2 ggplot aes .data geom_rect geom_hline labs
#'   facet_wrap vars scale_x_continuous scale_y_continuous
.plot_turnover <- function(df, time.col, type, panels, levels, ...) {
    p <- ggplot(df, aes(xmin = .data$.from, xmax = .data$.time,
            ymin = .data$.ymin, ymax = .data$.ymax, fill = .data$.part)) +
        # A thin outline separates adjacent intervals.
        do.call(geom_rect, modifyList(
            list(colour = "white", linewidth = 0.2), list(...))) +
        .scale_categorical("fill") +
        labs(x = time.col, fill = NULL,
            y = if (type == "gainloss") "features" else "turnover") +
        .theme_tempoPlot()
    if (type == "gainloss") {
        # Losses are drawn downwards but counted as positive numbers.
        p <- p + geom_hline(yintercept = 0, colour = "grey30") +
            scale_y_continuous(labels = abs)
    }
    if (!is.null(levels)) {
        p <- p + scale_x_continuous(breaks = seq_along(levels),
            labels = levels)
    }
    if (panels) {
        p <- p + facet_wrap(vars(.data$.panel))
    }
    p
}
