library(SummarizedExperiment)
library(SingleCellExperiment)

data("hitchip1006", package = "miaTime")

# Multi-sample subjects for order / gap tests.
.multi_df <- function(x = hitchip1006, n_subjects = 3L) {
    counts <- table(x$subject)
    ids <- names(counts)[counts >= 3L][seq_len(n_subjects)]
    keep <- x$subject %in% ids
    data.frame(
        subject = as.character(x$subject[keep]),
        time = as.numeric(x$time[keep]),
        sample = colnames(x)[keep],
        stringsAsFactors = FALSE
    )
}

test_that(".check_longitudinal accepts hitchip1006 and warns", {
    expect_warning(
        type <- .check_longitudinal(hitchip1006, "time", group = "subject"),
        "'group'"
    )
    expect_identical(type, "numeric")
})

test_that(".check_longitudinal validates class and columns", {
    expect_error(
        .check_longitudinal(list(), "time"),
        "'x'"
    )
    expect_error(
        .check_longitudinal(hitchip1006, "missing"),
        "'time.col'"
    )
    expect_error(
        .check_longitudinal(hitchip1006, "time", group = "missing"),
        "'group'"
    )
})

test_that(".check_longitudinal errors on NA and duplicates", {
    x_na <- hitchip1006
    colData(x_na)$time[1L] <- NA_real_
    expect_error(
        .check_longitudinal(x_na, "time", group = "subject"),
        "'time.col'"
    )
    x_na <- hitchip1006
    colData(x_na)$subject[1L] <- NA
    expect_error(
        .check_longitudinal(x_na, "time", group = "subject"),
        "'group' must not contain NA"
    )

    x_dup <- hitchip1006
    # Force a duplicated subject-time pair on a multi-sample subject.
    counts <- table(x_dup$subject)
    sid <- names(counts)[counts >= 2L][[1L]]
    i <- which(x_dup$subject == sid)[1:2]
    colData(x_dup)$time[i[[2L]]] <- colData(x_dup)$time[i[[1L]]]
    expect_error(
        suppressWarnings(
            .check_longitudinal(x_dup, "time", group = "subject")
        ),
        "duplicated"
    )

    # Without group: duplicated times across the object.
    x_time <- hitchip1006[, 1:5]
    colData(x_time)$time <- c(1, 1, 2, 3, 4)
    expect_error(
        .check_longitudinal(x_time, "time"),
        "duplicated"
    )
})

test_that(".time_type handles numeric, Date, factor and character", {
    expect_identical(.time_type(hitchip1006$time), "numeric")

    x_date <- hitchip1006
    colData(x_date)$time <- as.Date("2000-01-01") +
        as.integer(colData(x_date)$time)
    expect_identical(.time_type(x_date$time), "date")

    x_visit <- hitchip1006
    t <- colData(x_visit)$time
    colData(x_visit)$time <- factor(
        paste0("V", t),
        levels = paste0("V", sort(unique(t))),
        ordered = TRUE
    )
    expect_identical(.time_type(x_visit$time), "visit")

    x_fac <- hitchip1006
    colData(x_fac)$time <- factor(paste0("V", colData(x_fac)$time))
    expect_identical(.time_type(x_fac$time), "visit")

    x_chr <- hitchip1006
    colData(x_chr)$time <- as.character(colData(x_chr)$time)
    expect_error(.time_type(x_chr$time), "POSIXct or factor")
})

test_that("POSIXct times are continuous time in hours", {
    x <- hitchip1006[, 1:4]
    x$time <- as.POSIXct("2020-01-01", tz = "UTC") + c(0, 6, 30, 31) * 3600
    expect_identical(.check_longitudinal(x, "time"), "datetime")
    expect_equal(diff(.time_numeric(x$time)), c(6, 24, 1))

    df <- data.frame(subject = "A", time = x$time)
    expect_identical(.break_gaps(df, "time", "subject", max.gap = 12)$.segment,
        c(1L, 1L, 2L, 2L))

    # .time_restore() inverts .time_numeric(), keeping the time zone.
    tz <- as.POSIXct("2020-06-01 08:30", tz = "Europe/Helsinki")
    expect_identical(.time_restore(.time_numeric(tz), tz), tz)
    day <- as.Date("2020-06-01")
    expect_identical(.time_restore(.time_numeric(day), day), day)
    expect_identical(.time_restore(3, 1:2), 3)

    ev <- .check_events(data.frame(start = x$time[2L]), x, "time")
    expect_identical(ev$start, x$time[2L])
    expect_error(.check_events(data.frame(start = 6), x, "time"),
        "type of 'time.col'")
    expect_error(
        .check_events(data.frame(start = as.Date("2020-01-01")), x, "time"),
        "type of 'time.col'")
})

test_that(".order_time uses group and time from hitchip1006", {
    df <- .multi_df()
    # Shuffle rows so ordering is non-trivial.
    set.seed(1)
    df <- df[sample(nrow(df)), , drop = FALSE]
    out <- .order_time(df, "time", group = "subject")
    expect_identical(out$subject, sort(out$subject))
    # Within each subject, time is non-decreasing.
    split_out <- split(out$time, out$subject)
    expect_true(all(vapply(split_out, function(t) {
        all(diff(t) >= 0)
    }, logical(1))))

    # Without group: global time order.
    out2 <- .order_time(df, "time")
    expect_equal(out2$time, sort(df$time))

    # group must be used: subject-then-time differs from time-only.
    expect_false(identical(out$sample, out2$sample))
})

test_that(".break_gaps breaks numeric, Date and visit gaps", {
    df <- .multi_df(n_subjects = 2L)
    df <- .order_time(df, "time", group = "subject")

    # Break wherever consecutive times are more than 1 unit apart.
    gaps <- unlist(lapply(split(df$time, df$subject), diff), use.names = FALSE)
    expect_true(any(gaps > 1))
    out <- .break_gaps(df, "time", group = "subject", max.gap = 1)
    expect_true(max(out$.segment) >= 2L)

    out_null <- .break_gaps(
        df, "time", group = "subject", max.gap = NULL
    )
    expect_identical(out_null$.segment, rep(1L, nrow(df)))

    df_date <- df
    df_date$time <- as.Date("2000-01-01") + as.integer(df$time)
    out_date <- .break_gaps(
        df_date, "time", group = "subject", max.gap = 1
    )
    expect_identical(out_date$.segment, out$.segment)

    # Visit distance uses level steps, not numeric labels.
    df_visit <- data.frame(
        subject = c("A", "A", "A"),
        time = factor(
            c("V1", "V2", "V5"),
            levels = paste0("V", 1:5)
        ),
        stringsAsFactors = FALSE
    )
    out_visit <- .break_gaps(
        df_visit, "time", group = "subject", max.gap = 2
    )
    # V1->V2 = 1; V2->V5 = 3 > 2
    expect_identical(out_visit$.segment, c(1L, 1L, 2L))
    # Without max.gap, visits break at every skipped level, also when the
    # codes are numbers and visit says so.
    expect_identical(.break_gaps(df_visit, "time", "subject")$.segment,
        c(1L, 1L, 2L))
    df_codes <- transform(df_visit, time = as.integer(time))
    expect_identical(.break_gaps(df_codes, "time", "subject")$.segment,
        rep(1L, 3L))
    expect_identical(.break_gaps(df_codes, "time", "subject",
        visit = TRUE)$.segment, c(1L, 1L, 2L))

    expect_error(
        .break_gaps(df, "time", max.gap = -1),
        "'max.gap'"
    )
    expect_error(
        .break_gaps(df, "time", max.gap = c(1, 2)),
        "'max.gap'"
    )
})

test_that(".resolve_y finds diversity, taxon and PCA components", {
    x <- hitchip1006
    # Small subset keeps prcomp fast and stable.
    x <- x[, x$subject %in% names(table(x$subject))[1:20]]
    reducedDim(x, "PCA") <- prcomp(
        t(assay(x, "counts")),
        rank. = 2
    )$x[, 1:2, drop = FALSE]
    colnames(reducedDim(x, "PCA")) <- c("PC1", "PC2")

    y_cd <- .resolve_y(x, "diversity")
    expect_equal(as.numeric(y_cd), as.numeric(x$diversity))
    expect_identical(names(y_cd), colnames(x))

    taxon <- rownames(x)[[1L]]
    y_feat <- .resolve_y(x, taxon, assay.type = "counts")
    expect_equal(as.numeric(y_feat), as.numeric(assay(x, "counts")[taxon, ]))

    y_idx <- .resolve_y(x, "PCA:1")
    expect_equal(
        as.numeric(y_idx),
        as.numeric(reducedDim(x, "PCA")[, 1L])
    )

    y_name <- .resolve_y(x, "PCA:PC2")
    expect_equal(
        as.numeric(y_name),
        as.numeric(reducedDim(x, "PCA")[, "PC2"])
    )
})

test_that(".resolve_y errors on ambiguous, missing and non-numeric", {
    x <- hitchip1006[, 1:30]
    taxon <- rownames(x)[[1L]]
    colData(x)[[taxon]] <- seq_len(ncol(x))
    expect_error(.resolve_y(x, taxon), "more than one place")

    expect_error(.resolve_y(x, "missing_y"), "not found")

    expect_error(.resolve_y(x, "sex"), "numeric")
    expect_error(.resolve_y(hitchip1006, taxon, assay.type = "clr"),
        "'assay.type' must be one of assayNames")
})
