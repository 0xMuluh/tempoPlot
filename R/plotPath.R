#' Plot trajectories of subjects in two dimensions
#'
#' \code{plotPath()} joins each subject's samples in time order in a
#' two-dimensional space: two components of an ordination, or two variables
#' such as features or colData columns.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date} or
#'   a factor of visits.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param dimred \code{Character scalar}. Specifies the
#'   \code{reducedDim(x)} to plot. (Default: \code{NULL})
#'
#' @param ncomponents \code{Numeric} or \code{Character vector}. The two
#'   components of \code{dimred} to plot, by index or name.
#'   (Default: \code{c(1, 2)})
#'
#' @param x.var,y.var \code{Character scalar}. Variables to plot instead of
#'   \code{dimred}: a column of \code{colData(x)} or a feature of
#'   \code{assay.type}. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay that holds
#'   features named in \code{x.var} or \code{y.var}. (Default: \code{NULL})
#'
#' @param origin \code{Character scalar}. Whether paths start where the
#'   first sample lies (\code{"none"}) or at zero, by subtracting each
#'   subject's first sample (\code{"baseline"}). (Default: \code{"none"})
#'
#' @param arrows \code{Logical scalar}. Whether paths end in an arrow
#'   pointing forward in time. (Default: \code{TRUE})
#'
#' @param colour.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that colours the paths. (Default: \code{NULL})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that groups subjects for \code{centroids}.
#'   (Default: \code{NULL})
#'
#' @param centroids \code{Logical scalar}. Draws the mean path of each
#'   \code{group.by} group, or of all subjects: the mean position at each
#'   time. Best with times shared by the subjects, such as visits.
#'   (Default: \code{FALSE})
#'
#' @param field \code{Integer scalar}. Draws the mean movement to the next
#'   sample from each cell of a \code{field} by \code{field} grid as an
#'   arrow. If \code{NULL}, draws no field. (Default: \code{NULL})
#'
#' @param lag \code{Integer scalar}. Instead of paths, draws the value of
#'   \code{x.var}, or the first component of \code{dimred}, against its
#'   value \code{lag} samples later in the same subject, with the identity
#'   line. (Default: \code{NULL})
#'
#' @param animate \code{Logical scalar}. Reveals the paths over time with
#'   \code{\link[gganimate:transition_reveal]{gganimate::transition_reveal()}};
#'   requires gganimate. (Default: \code{FALSE})
#'
#' @param max.gap \code{Numeric scalar}. Breaks a path where consecutive
#'   samples of a subject are further apart than this, in \code{time.col}
#'   units. \code{NULL} never breaks paths. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_path]{geom_path}}.
#'
#' @details
#' Give either \code{dimred} or both \code{x.var} and \code{y.var}; with
#' \code{lag}, \code{x.var} alone is enough. Points near the identity line
#' of a lag plot change little over \code{lag} samples.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @name plotPath
#'
#' @examples
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
#' # Subjects sampled five times
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]
#'
#' # Paths of subjects in a PCoA of Bray-Curtis dissimilarities
#' tse <- mia::transformAssay(tse, method = "relabundance")
#' tse <- mia::addMDS(tse, assay.type = "relabundance", method = "bray",
#'     name = "PCoA")
#' plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
#'     colour.by = "subject")
#'
#' # The same paths, all starting at zero, to compare their directions
#' plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
#'     origin = "baseline", colour.by = "subject")
#'
#' # Two genera against each other
#' plotPath(tse, time.col = "time", group = "subject",
#'     x.var = "Bacteroides vulgatus et rel.",
#'     y.var = "Prevotella melaninogenica et rel.",
#'     assay.type = "relabundance", colour.by = "subject")
#'
#' # Mean movement over a 6 x 6 grid of the PCoA
#' plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
#'     field = 6, arrows = FALSE)
#'
#' # Mean path per visit and sex
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' plotPath(tse, time.col = "visit", group = "subject", dimred = "PCoA",
#'     group.by = "sex", centroids = TRUE)
#'
#' # Diversity against diversity at the next sample
#' plotPath(tse, time.col = "time", group = "subject", x.var = "diversity",
#'     lag = 1)
NULL

#' @rdname plotPath
#' @importFrom utils modifyList
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom methods setMethod
#' @export
setMethod("plotPath", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, dimred = NULL, ncomponents = c(1, 2),
            x.var = NULL, y.var = NULL, assay.type = NULL,
            origin = c("none", "baseline"), arrows = TRUE, colour.by = NULL,
            group.by = NULL, centroids = FALSE, field = NULL, lag = NULL,
            animate = FALSE, max.gap = NULL, ...) {
        origin <- match.arg(origin)
        for (arg in c("arrows", "centroids", "animate")) {
            if (!(isTRUE(get(arg)) || isFALSE(get(arg)))) {
                stop("'", arg, "' must be TRUE or FALSE.", call. = FALSE)
            }
        }
        for (arg in c("field", "lag")) {
            v <- get(arg)
            if (!is.null(v) && !(is.numeric(v) && length(v) == 1L &&
                    v >= 1)) {
                stop("'", arg, "' must be NULL or a single number of at ",
                    "least 1.", call. = FALSE)
            }
        }
        if (!is.null(lag) && (centroids || !is.null(field) || animate ||
                origin != "none")) {
            stop("'lag' cannot be combined with 'centroids', 'field', ",
                "'animate' or 'origin'.", call. = FALSE)
        }
        if (animate && !requireNamespace("gganimate", quietly = TRUE)) {
            stop("'animate' requires the gganimate package.", call. = FALSE)
        }
        .check_longitudinal(x, time.col, group)
        # Errors on character and other unsupported time types.
        .time_type(x, time.col)
        .check_column(x, colour.by)
        .check_column(x, group.by)
        if (!is.null(lag) && is.null(dimred) && !is.null(x.var)) {
            y.var <- x.var
        }
        if (is.null(dimred) == (is.null(x.var) || is.null(y.var))) {
            stop("Give either 'dimred' or both 'x.var' and 'y.var'.",
                call. = FALSE)
        }

        cols <- unique(c(time.col, group, colour.by, group.by))
        df <- as.data.frame(colData(x)[, cols, drop = FALSE],
            optional = TRUE)
        if (is.null(dimred)) {
            df$.x <- unname(.resolve_y(x, x.var, assay.type, "x.var"))
            df$.y <- unname(.resolve_y(x, y.var, assay.type, "y.var"))
            labels <- c(x.var, y.var)
        } else {
            labels <- paste0(dimred, ":", ncomponents)
            hits <- lapply(labels, .match_dimred_component, x = x)
            if (length(labels) != 2L || any(vapply(hits, is.null,
                    logical(1)))) {
                stop("'dimred' must be in reducedDimNames(x) and ",
                    "'ncomponents' must give two of its components.",
                    call. = FALSE)
            }
            rd <- reducedDim(x, dimred)
            df$.x <- as.numeric(rd[, hits[[1L]]$j])
            df$.y <- as.numeric(rd[, hits[[2L]]$j])
        }
        df <- .order_time(df, time.col, group)
        df <- .break_gaps(df, time.col, group, max.gap)
        if (origin == "baseline") {
            first <- function(v) ave(v, df[[group]], FUN = function(u) u[1L])
            df$.x <- df$.x - first(df$.x)
            df$.y <- df$.y - first(df$.y)
        }
        df$.line <- interaction(df[[group]], df$.segment, drop = TRUE)
        df$.t <- .time_numeric(df[[time.col]])

        if (!is.null(lag)) {
            return(.plot_lag(df, group, labels[1L], as.integer(lag),
                colour.by, ...))
        }
        p <- .plot_path(df, labels, arrows, colour.by,
            fade = centroids || !is.null(field), ...)
        if (!is.null(field)) {
            p <- p + .field_layer(df, as.integer(field))
        }
        if (centroids) {
            p <- p + .centroid_layer(df, time.col, group.by, arrows)
        }
        if (animate) {
            p <- p + gganimate::transition_reveal(.data$.t)
        }
        p
    }
)

#' @importFrom ggplot2 ggplot aes .data geom_path geom_point labs arrow unit
.plot_path <- function(df, labels, arrows, colour.by, fade = FALSE, ...) {
    alpha <- if (fade) list(alpha = 0.3) else list()
    ggplot(df, aes(x = .data$.x, y = .data$.y,
            colour = if (!is.null(colour.by)) .data[[colour.by]])) +
        do.call(geom_path, c(list(aes(group = .data$.line),
            arrow = if (arrows) arrow(length = unit(0.2, "cm"))),
            modifyList(alpha, list(...)))) +
        do.call(geom_point, c(list(size = 1), alpha)) +
        labs(x = labels[1L], y = labels[2L], colour = colour.by)
}

# Mean path of each group.by group (or of all subjects) over time.
#' @importFrom ggplot2 geom_path aes .data arrow unit
#' @importFrom stats aggregate
.centroid_layer <- function(df, time.col, group.by, arrows) {
    by <- if (is.null(group.by)) rep("", nrow(df)) else
        as.character(df[[group.by]])
    cen <- aggregate(cbind(.x, .y) ~ .by + .t,
        data = data.frame(.x = df$.x, .y = df$.y, .by = by,
            .t = .time_numeric(df[[time.col]])), FUN = mean)
    cen <- cen[order(cen$.by, cen$.t), , drop = FALSE]
    if (is.null(group.by)) {
        return(geom_path(aes(x = .data$.x, y = .data$.y, group = .data$.by),
            data = cen, inherit.aes = FALSE, linewidth = 1.2,
            colour = "black",
            arrow = if (arrows) arrow(length = unit(0.3, "cm"))))
    }
    list(geom_path(aes(x = .data$.x, y = .data$.y, group = .data$.by,
            linetype = .data$.by), data = cen, inherit.aes = FALSE,
            linewidth = 1.2, colour = "black",
            arrow = if (arrows) arrow(length = unit(0.3, "cm"))),
        labs(linetype = group.by))
}

# Mean movement to the next sample of the same path segment, from each
# cell of a bins x bins grid, drawn from the cell centre.
#' @importFrom ggplot2 geom_segment aes .data arrow unit
#' @importFrom stats aggregate
.field_layer <- function(df, bins) {
    n <- nrow(df)
    i <- which(df$.line[-1L] == df$.line[-n])
    if (length(i) == 0L) {
        return(list())
    }
    cell <- function(v) {
        br <- seq(min(v), max(v), length.out = bins + 1L)
        mid <- (br[-1L] + br[-length(br)]) / 2
        mid[findInterval(v, br, rightmost.closed = TRUE, all.inside = TRUE)]
    }
    move <- data.frame(.cx = cell(df$.x)[i], .cy = cell(df$.y)[i],
        .dx = df$.x[i + 1L] - df$.x[i], .dy = df$.y[i + 1L] - df$.y[i])
    move <- aggregate(cbind(.dx, .dy) ~ .cx + .cy, data = move, FUN = mean)
    geom_segment(aes(x = .data$.cx, y = .data$.cy,
            xend = .data$.cx + .data$.dx, yend = .data$.cy + .data$.dy),
        data = move, inherit.aes = FALSE, colour = "grey20",
        arrow = arrow(length = unit(0.15, "cm")))
}

# Value at each sample against the value lag samples later in the same
# subject, with the identity line.
#' @importFrom ggplot2 ggplot aes .data geom_point geom_abline labs
.plot_lag <- function(df, group, label, lag, colour.by, ...) {
    n <- nrow(df)
    if (n <= lag) {
        stop("No subject has more than 'lag' samples.", call. = FALSE)
    }
    i <- seq_len(n - lag)
    i <- i[df[[group]][i] == df[[group]][i + lag]]
    if (length(i) == 0L) {
        stop("No subject has more than 'lag' samples.", call. = FALSE)
    }
    pairs <- df[i, , drop = FALSE]
    pairs$.later <- df$.x[i + lag]
    ggplot(pairs, aes(x = .data$.x, y = .data$.later,
            colour = if (!is.null(colour.by)) .data[[colour.by]])) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
            colour = "grey40") +
        geom_point(...) +
        labs(x = paste(label, "at a sample"),
            y = paste(label, lag, if (lag == 1L) "sample" else "samples",
                "later"),
            colour = colour.by)
}
