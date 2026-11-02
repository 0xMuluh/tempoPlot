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
    cd_names <- names(colData(x))
    if (!.is_one_name(time.col) || !time.col %in% cd_names) {
        stop("'time.col' must be one column name of colData(x).",
            call. = FALSE)
    }
    if (!is.null(group) &&
            (!.is_one_name(group) || !group %in% cd_names)) {
        stop("'group' must be NULL or one column name of colData(x).",
            call. = FALSE)
    }

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
# Places checked: colData column, assay row, "<dimred>:<component>".
#' @importFrom SummarizedExperiment assay colData
.resolve_y <- function(x, y, assay.type = NULL) {
    if (!.is_one_name(y)) {
        stop("'y' must be a single column, feature or dimred name.",
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
            "'y' was found in more than one place: ",
            paste(places, collapse = ", "),
            ".",
            call. = FALSE
        )
    }
    if (length(places) == 0L) {
        stop("'y' was not found in colData, assay rows or reducedDim ",
            "components.", call. = FALSE)
    }

    if (identical(places, "colData")) {
        col <- colData(x)[[y]]
        if (!is.numeric(col)) {
            stop("'y' must be a numeric column of colData(x).",
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
