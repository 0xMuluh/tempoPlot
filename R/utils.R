# Internal helpers for longitudinal plot functions.
# Not exported; documented here with comments only.

# Check that x is longitudinal-ready for plotting.
# Errors on bad class, missing columns, NA times and duplicated
# subject-time pairs; warns on single-sample subjects when group is set.
# Returns x invisibly.
#' @importFrom SummarizedExperiment colData
.check_longitudinal <- function(x, time.col, group = NULL) {
    if (!inherits(x, "SummarizedExperiment")) {
        stop("'x' must be a SummarizedExperiment.", call. = FALSE)
    }
    if (!.is_one_name(time.col) || !time.col %in% names(colData(x))) {
        stop("'time.col' must be one column name of colData(x).",
            call. = FALSE)
    }
    .check_column(x, group)

    time <- colData(x)[[time.col]]
    if (anyNA(time)) {
        stop("'time.col' must not contain NA values.", call. = FALSE)
    }

    if (is.null(group)) {
        if (anyDuplicated(time)) {
            stop("'time.col' must not contain duplicated times.",
                call. = FALSE)
        }
    } else {
        subject <- colData(x)[[group]]
        if (anyDuplicated(data.frame(subject, time))) {
            stop("'group' and 'time.col' must not contain duplicated ",
                "subject-time pairs.", call. = FALSE)
        }
        counts <- table(subject)
        singles <- names(counts)[counts == 1L]
        if (length(singles) > 0L) {
            warning(
                "'group' has subjects with a single sample: ",
                paste(singles, collapse = ", "),
                ".",
                call. = FALSE
            )
        }
    }
    invisible(x)
}

# Classify the time column as "numeric", "date" or "visit".
.time_type <- function(x, time.col) {
    time <- colData(x)[[time.col]]
    if (is.numeric(time)) {
        return("numeric")
    }
    if (inherits(time, "Date")) {
        return("date")
    }
    if (is.factor(time)) {
        return("visit")
    }
    stop(
        "'time.col' must be numeric, Date or factor; ",
        "convert character times to a factor.",
        call. = FALSE
    )
}

# Order a data.frame by group (optional) then time.
# Factor times follow level order; ties keep original order.
.order_time <- function(df, time.col, group = NULL) {
    if (is.null(group)) {
        ord <- order(df[[time.col]])
    } else {
        ord <- order(df[[group]], df[[time.col]])
    }
    df[ord, , drop = FALSE]
}

# Time as numbers: Date as days, factor as level codes. The type is
# validated by .time_type().
.time_numeric <- function(time) {
    if (is.factor(time)) as.integer(time) else as.numeric(time)
}

# Add integer column .segment, increasing within each subject when the
# gap to the previous sample exceeds max.gap. Distance is numeric
# difference, days for Date, or level steps for visit factors.
# max.gap = NULL keeps one segment per subject.
#' @importFrom stats ave
.break_gaps <- function(df, time.col, group = NULL, max.gap = NULL) {
    if (!is.null(max.gap) &&
            !(is.numeric(max.gap) && length(max.gap) == 1L &&
                is.finite(max.gap) && max.gap > 0)) {
        stop("'max.gap' must be NULL or a single positive number.",
            call. = FALSE)
    }

    n <- nrow(df)
    if (n == 0L || is.null(max.gap)) {
        df$.segment <- rep(1L, n)
        return(df)
    }

    if (is.null(group)) {
        ids <- rep(1L, n)
    } else {
        ids <- df[[group]]
    }
    tnum <- .time_numeric(df[[time.col]])
    new_subject <- c(TRUE, ids[-1L] != ids[-n])
    large_gap <- c(FALSE, diff(tnum) > max.gap)
    flag <- as.integer(new_subject | large_gap)
    df$.segment <- as.integer(ave(flag, ids, FUN = cumsum))
    df
}

# Resolve y to a numeric vector of length ncol(x), named by colnames(x).
# Places checked: colData column, assay row, "<dimred>:<component>". Error
# messages name the argument arg.
#' @importFrom SummarizedExperiment assay colData
.resolve_y <- function(x, y, assay.type = NULL, arg = "y") {
    if (!.is_one_name(y)) {
        stop("'", arg, "' must be a single column, feature or dimred name.",
            call. = FALSE)
    }

    places <- character(0)
    if (y %in% names(colData(x))) {
        places <- c(places, "colData")
    }
    if (y %in% rownames(x)) {
        places <- c(places, "assay")
    }
    dimred_match <- .match_dimred_component(x, y)
    if (!is.null(dimred_match)) {
        places <- c(places, "reducedDim")
    }

    if (length(places) > 1L) {
        stop(
            "'", arg, "' was found in more than one place: ",
            paste(places, collapse = ", "),
            ".",
            call. = FALSE
        )
    }
    if (length(places) == 0L) {
        stop("'", arg, "' was not found in colData, assay rows or reducedDim ",
            "components.", call. = FALSE)
    }

    if (identical(places, "colData")) {
        col <- colData(x)[[y]]
        if (!is.numeric(col)) {
            stop("'", arg, "' must be a numeric column of colData(x).",
                call. = FALSE)
        }
        value <- as.numeric(col)
    } else if (identical(places, "assay")) {
        value <- as.numeric(assay(x, assay.type %||% 1L)[y, ])
    } else {
        rd <- reducedDim(x, dimred_match$dimred)
        value <- as.numeric(rd[, dimred_match$j])
    }

    names(value) <- colnames(x)
    value
}

# TRUE when x is a single non-NA character string.
.is_one_name <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x)
}

# Error unless col is NULL or one column name of colData(x) or rowData(x).
# The message names the argument passed as col.
#' @importFrom SummarizedExperiment colData rowData
.check_column <- function(x, col, arg = deparse(substitute(col)),
        where = c("colData", "rowData")) {
    where <- match.arg(where)
    cols <- names(if (where == "colData") colData(x) else rowData(x))
    if (!is.null(col) && !(.is_one_name(col) && col %in% cols)) {
        stop("'", arg, "' must be NULL or one column name of ", where,
            "(x).", call. = FALSE)
    }
}

# If y is "<dimred>:<component>", return list(dimred, j) for the matching
# column; otherwise NULL. Split at the first ':' only. Component may be
# a column name or a 1-based index.
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
.match_dimred_component <- function(x, y) {
    if (!inherits(x, "SingleCellExperiment")) {
        return(NULL)
    }
    i <- regexpr(":", y, fixed = TRUE)
    if (i < 1L) {
        return(NULL)
    }
    dimred <- substring(y, 1L, i - 1L)
    component <- substring(y, i + 1L)
    if (!nzchar(dimred) || !nzchar(component) ||
            !dimred %in% reducedDimNames(x)) {
        return(NULL)
    }
    rd <- reducedDim(x, dimred)
    if (component %in% colnames(rd)) {
        return(list(dimred = dimred, j = component))
    }
    if (grepl("^[0-9]+$", component)) {
        idx <- as.integer(component)
        if (idx >= 1L && idx <= ncol(rd)) {
            return(list(dimred = dimred, j = idx))
        }
    }
    NULL
}

# Check events and return them as a data.frame, or NULL. Columns: start
# (required), end and label (optional) and group (optional; without it,
# each event applies to all subjects). start and end must have the type
# of time.col: numeric, Date, or levels of the visit factor (returned as
# that factor). end NA marks a single time point.
.check_events <- function(events, x, time.col, group = NULL) {
    if (is.null(events)) {
        return(NULL)
    }
    if (!is.data.frame(events) || !"start" %in% names(events)) {
        stop("'events' must be a data.frame with a column 'start'.",
            call. = FALSE)
    }
    events <- as.data.frame(events)
    time <- colData(x)[[time.col]]
    for (col in intersect(c("start", "end"), names(events))) {
        v <- events[[col]]
        if (is.factor(time)) {
            ok <- all(is.na(v) | as.character(v) %in% levels(time))
            v <- factor(as.character(v), levels = levels(time))
        } else if (inherits(time, "Date")) {
            ok <- inherits(v, "Date")
        } else {
            ok <- is.numeric(v)
        }
        if (!ok) {
            stop("'events$", col, "' must have the type of 'time.col'",
                if (is.factor(time)) " and contain its levels", ".",
                call. = FALSE)
        }
        events[[col]] <- v
    }
    if (anyNA(events$start)) {
        stop("'events$start' must not contain NA values.", call. = FALSE)
    }
    if ("end" %in% names(events) && any(
            .time_numeric(events$end) < .time_numeric(events$start),
            na.rm = TRUE)) {
        stop("'events$end' must not be before 'events$start'.",
            call. = FALSE)
    }
    if ("label" %in% names(events)) {
        events$label <- as.character(events$label)
    }
    events
}

# Layers for events checked by .check_events(): bands (rows with end) or
# markers (rows without end). Events with a group column are drawn on the
# row of their subject when y.levels gives the subject rows, bottom to
# top; the others span the whole plot. label sets the fill, or the
# linetype of markers that span the plot. Returns a
# list of layers, empty without events of that kind.
#' @importFrom ggplot2 aes .data geom_point geom_rect geom_vline
.events_layer <- function(events, what = c("bands", "markers"),
        group = NULL, y.levels = NULL) {
    what <- match.arg(what)
    if (is.null(events)) {
        return(list())
    }
    band <- if ("end" %in% names(events)) !is.na(events$end) else FALSE
    keep <- if (what == "bands") band else !band
    ev <- events[rep_len(keep, nrow(events)), , drop = FALSE]
    rows <- !is.null(y.levels) && !is.null(group) && group %in% names(ev)
    if (rows) {
        ev$.y <- match(as.character(ev[[group]]), y.levels)
        ev <- ev[!is.na(ev$.y), , drop = FALSE]
    }
    if (nrow(ev) == 0L) {
        return(list())
    }
    # Visits are drawn at their level codes.
    pos <- function(v) if (is.factor(v)) as.integer(v) else v
    ev$.start <- pos(ev$start)
    fill <- if ("label" %in% names(ev)) "label" else NULL
    if (what == "bands") {
        ev$.end <- pos(ev$end)
        ev$.ymin <- if (rows) ev$.y - 0.4 else -Inf
        ev$.ymax <- if (rows) ev$.y + 0.4 else Inf
        mapping <- aes(xmin = .data$.start, xmax = .data$.end,
            ymin = .data$.ymin, ymax = .data$.ymax)
        if (is.null(fill)) {
            return(list(geom_rect(mapping, data = ev, inherit.aes = FALSE,
                fill = "grey50", alpha = 0.3)))
        }
        mapping <- aes(xmin = .data$.start, xmax = .data$.end,
            ymin = .data$.ymin, ymax = .data$.ymax, fill = .data$label)
        return(list(geom_rect(mapping, data = ev, inherit.aes = FALSE,
            alpha = 0.3)))
    }
    if (!rows) {
        if (is.null(fill)) {
            return(list(geom_vline(aes(xintercept = .data$.start),
                data = ev, linetype = "dashed", colour = "grey30")))
        }
        return(list(geom_vline(
            aes(xintercept = .data$.start, linetype = .data$label),
            data = ev, colour = "grey30")))
    }
    ev$.y <- ev$.y + 0.3
    if (is.null(fill)) {
        return(list(geom_point(aes(x = .data$.start, y = .data$.y),
            data = ev, inherit.aes = FALSE, shape = 25, size = 2.5,
            fill = "grey30")))
    }
    list(geom_point(aes(x = .data$.start, y = .data$.y, fill = .data$label),
        data = ev, inherit.aes = FALSE, shape = 25, size = 2.5))
}

# Sample-by-sample dissimilarity matrix stored by mia::addDissimilarity()
# in metadata(x)[[dissimilarity]], as a full matrix in colnames(x) order.
#' @importFrom S4Vectors metadata
.get_dissimilarity <- function(x, dissimilarity) {
    if (!.is_one_name(dissimilarity) ||
            !dissimilarity %in% names(metadata(x))) {
        stop("'dissimilarity' must name a matrix in metadata(x), as stored ",
            "by mia::addDissimilarity().", call. = FALSE)
    }
    d <- as.matrix(metadata(x)[[dissimilarity]])
    if (!is.numeric(d) || nrow(d) != ncol(x) || ncol(d) != ncol(x)) {
        stop("metadata(x)$", dissimilarity, " must be a sample by sample ",
            "dissimilarity.", call. = FALSE)
    }
    if (!is.null(rownames(d)) && !is.null(colnames(x))) {
        if (!setequal(rownames(d), colnames(x))) {
            stop("metadata(x)$", dissimilarity, " must be named by ",
                "colnames(x).", call. = FALSE)
        }
        d <- d[colnames(x), colnames(x)]
    }
    d
}
