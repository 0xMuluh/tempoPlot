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
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects. If \code{NULL}, all samples
#'   form one series. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay; a
#'   feature is present in a sample when its value is above zero.
#'   (Default: \code{"counts"})
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
#' plotTurnover(minimalgut, time.col = "Time.hr", group = "StudyIdentifier")
#'
#' # Turnover components, averaged over the bioreactors
#' plotTurnover(minimalgut, time.col = "Time.hr", group = "StudyIdentifier",
#'     type = "components", level = "group")
NULL

#' @rdname plotTurnover
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom stats aggregate ave
#' @importFrom methods setMethod
#' @export
setMethod("plotTurnover", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group = NULL, assay.type = "counts",
            type = c("gainloss", "components"),
            level = c("subject", "group"), group.by = NULL, ...) {
        type <- match.arg(type)
        level <- match.arg(level)
        .check_longitudinal(x, time.col, group)
        ttype <- .time_type(x, time.col)
        .check_column(x, group.by)
        if (!.is_one_name(assay.type) || !assay.type %in% assayNames(x)) {
            stop("'assay.type' must be one of assayNames(x).",
                call. = FALSE)
        }

        cols <- unique(c(time.col, group, group.by))
        cd <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        cd$.i <- seq_len(ncol(x))
        cd <- .order_time(cd, time.col, group)
        n <- nrow(cd)
        subject <- if (is.null(group)) rep("", n) else
            as.character(cd[[group]])
        same <- which(subject[-1L] == subject[-n])
        if (length(same) == 0L) {
            stop("No subject has more than one sample.", call. = FALSE)
        }
        present <- as.matrix(assay(x, assay.type)) > 0
        before <- present[, cd$.i[same], drop = FALSE]
        after <- present[, cd$.i[same + 1L], drop = FALSE]
        gained <- colSums(!before & after, na.rm = TRUE)
        lost <- colSums(before & !after, na.rm = TRUE)
        total <- colSums(before | after, na.rm = TRUE)

        panel <- if (level == "subject") subject[same + 1L] else if
            (is.null(group.by)) rep("", length(same)) else
            as.character(cd[[group.by]][same + 1L])
        tnum <- .time_numeric(cd[[time.col]])
        time <- tnum[same + 1L]
        from <- tnum[same]
        df <- if (type == "gainloss") {
            data.frame(.panel = panel, .from = from, .time = time,
                .part = rep(c("gained", "lost"), each = length(same)),
                .value = c(gained, -lost))
        } else {
            share <- function(v) ifelse(total > 0, v / total, 0)
            data.frame(.panel = panel, .from = from, .time = time,
                .part = rep(c("appearance", "disappearance"),
                    each = length(same)),
                .value = c(share(gained), share(lost)))
        }
        if (level == "group") {
            df <- aggregate(cbind(.from, .value) ~ .panel + .time + .part,
                data = df, FUN = mean)
        }
        # Gains and appearance from zero up; losses from zero down and
        # disappearance stacked on appearance.
        df <- df[order(df$.panel, df$.time, df$.part), , drop = FALSE]
        df$.ymin <- pmin(df$.value, 0)
        df$.ymax <- pmax(df$.value, 0)
        if (type == "components") {
            up <- ave(df$.value, df$.panel, df$.time, FUN = function(v)
                c(0, cumsum(v)[-length(v)]))
            df$.ymin <- up
            df$.ymax <- up + df$.value
        }
        if (ttype == "visit") {
            df$.from <- df$.from + 0.05
            df$.time <- df$.time - 0.05
        }
        if (ttype == "date") {
            df$.from <- as.Date(df$.from, origin = "1970-01-01")
            df$.time <- as.Date(df$.time, origin = "1970-01-01")
        }
        levels <- if (ttype == "visit") levels(cd[[time.col]])
        panels <- level == "group" && !is.null(group.by) ||
            level == "subject" && !is.null(group)
        .plot_turnover(df, time.col, type, panels, levels, ...)
    }
)

#' @importFrom utils modifyList
#' @importFrom ggplot2 ggplot aes .data geom_rect geom_hline labs
#'   facet_wrap vars scale_fill_manual scale_x_continuous
.plot_turnover <- function(df, time.col, type, panels, levels, ...) {
    colours <- c(gained = "#2166ac", lost = "#b2182b",
        appearance = "#2166ac", disappearance = "#b2182b")
    p <- ggplot(df, aes(xmin = .data$.from, xmax = .data$.time,
            ymin = .data$.ymin, ymax = .data$.ymax, fill = .data$.part)) +
        # A thin outline separates adjacent intervals.
        do.call(geom_rect, modifyList(
            list(colour = "white", linewidth = 0.2), list(...))) +
        scale_fill_manual(values = colours) +
        labs(x = time.col, fill = NULL,
            y = if (type == "gainloss") "features" else "turnover")
    if (type == "gainloss") {
        p <- p + geom_hline(yintercept = 0, colour = "grey30")
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
