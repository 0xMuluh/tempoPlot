library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

sce <- SingleCellExperiment(
    assays = list(counts = matrix(c(0, 1, 2, 3, 0, 5), 2L,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3")))),
    colData = DataFrame(id = c("a", "a", "b")),
    reducedDims = list(PCoA = matrix(1:9, 3L,
        dimnames = list(NULL, c("PC1", "PC2", "PC3")))))

test_that(".get_assay returns the named assay as a matrix", {
    expect_identical(.get_assay(sce, "counts"), assay(sce, "counts"))
    for (bad in list("none", NULL, 1L, c("counts", "counts"))) {
        expect_error(.get_assay(sce, bad),
            "'assay.type' must be one of assayNames\\(x\\)")
    }
})

test_that(".require_column accepts one colData column and names the arg", {
    expect_no_error(.require_column(sce, "id"))
    group <- NULL
    expect_error(.require_column(sce, group),
        "'group' must be one column name of colData\\(x\\)")
    expect_error(.require_column(sce, "none", "state.col"), "'state.col'")
    expect_error(.require_column(sce, c("id", "id"), "group"), "'group'")
})

test_that(".get_reduced_dim takes a count or a vector of components", {
    rd <- .get_reduced_dim(sce, "PCoA", 2)
    expect_identical(colnames(rd), c("PCoA:1", "PCoA:2"))
    expect_equal(unname(rd), cbind(1:3, 4:6))
    rd <- .get_reduced_dim(sce, "PCoA", c(3, 1), n = 2L)
    expect_equal(unname(rd), cbind(7:9, 1:3))
    rd <- .get_reduced_dim(sce, "PCoA", c("PC2", "PC3"))
    expect_equal(unname(rd), cbind(4:6, 7:9))
    expect_identical(dim(.get_reduced_dim(sce, "PCoA", 3)), c(3L, 3L))
})

test_that(".get_reduced_dim errors on wrong dimred or components", {
    msg <- "'dimred' must be one of reducedDimNames"
    expect_error(.get_reduced_dim(sce, "none", c(1, 2)), msg)
    expect_error(.get_reduced_dim(sce, NULL, c(1, 2)), msg)
    expect_error(.get_reduced_dim(as(sce, "SummarizedExperiment"), "PCoA",
        c(1, 2)), msg)
    expect_error(.get_reduced_dim(sce, "PCoA", c(1, 4)), "'ncomponents'")
    expect_error(.get_reduced_dim(sce, "PCoA", 3, n = 2L), "select 2")
    expect_error(.get_reduced_dim(sce, "none", c(1, 2), arg = "dimred2"),
        "'dimred2'")
})

test_that(".primary_names maps experiment columns to primary samples", {
    e <- function(cols) SummarizedExperiment(
        assays = list(counts = matrix(0, 1L, length(cols),
            dimnames = list("f1", cols))))
    mae <- MultiAssayExperiment(
        experiments = ExperimentList(a = e(c("a1", "a2", "a3")),
            b = e(c("b3", "b1"))),
        colData = DataFrame(row.names = c("p1", "p2", "p3")),
        sampleMap = DataFrame(
            assay = factor(c("a", "a", "a", "b", "b")),
            primary = c("p1", "p2", "p3", "p1", "p3"),
            colname = c("a1", "a2", "a3", "b1", "b3")))
    expect_identical(.primary_names(mae, "a"), c("p1", "p2", "p3"))
    expect_identical(.primary_names(mae, "b"), c("p3", "p1"))
})

test_that(".next_in_subject pairs rows of the same subject", {
    ids <- c("a", "a", "b", "b", "b")
    expect_identical(.next_in_subject(ids), c(1L, 3L, 4L))
    expect_identical(.next_in_subject(ids, lag = 2L), 3L)
    expect_identical(.next_in_subject(c("a", NA, NA, "b")), integer(0))
    expect_identical(.next_in_subject(c("a", "a", NA)), 1L)
    expect_identical(.next_in_subject("a"), integer(0))
    expect_identical(.next_in_subject(character(0)), integer(0))
    # Pairs across a gap segment are skipped.
    expect_identical(.next_in_subject(ids, segment = c(1, 1, 1, 2, 2)),
        c(1L, 4L))
})

test_that(".check_flag and .check_number name the argument and range", {
    flag <- NA
    expect_error(.check_flag(flag), "'flag' must be TRUE or FALSE.",
        fixed = TRUE)
    expect_silent(.check_flag(FALSE))
    n <- 0
    expect_silent(.check_number(n, 0))
    expect_error(.check_number(n, 1),
        "'n' must be a single number of at least 1.", fixed = TRUE)
    expect_error(.check_number(n, 0, open = TRUE, null = TRUE),
        "'n' must be NULL or a single number above 0.", fixed = TRUE)
    expect_silent(.check_number(NULL, null = TRUE))
    expect_error(.check_number(n, 0, 1, open = TRUE),
        "'n' must be a single number in (0, 1].", fixed = TRUE)
    expect_silent(.check_number(1, 0, 1, open = TRUE))
    for (bad in list("a", c(1, 2), NA_real_, NULL)) {
        expect_error(.check_number(bad, arg = "x"),
            "'x' must be a single number.", fixed = TRUE)
    }
})

test_that(".present counts values above detection", {
    m <- matrix(c(0, 1, 2, NA), 2L)
    expect_identical(.present(m), m > 0)
    expect_identical(.present(m, 1), matrix(c(FALSE, FALSE, TRUE, NA), 2L))
    for (bad in list("a", c(1, 2), NA_real_, NULL)) {
        expect_error(.present(m, bad), "'detection' must be a single number")
    }
})

test_that(".tile_edges reaches halfway to neighbours, visits +-0.5", {
    expect_identical(.tile_edges(c(0, 1, 3)),
        list(lower = c(-0.5, 0.5, 2), upper = c(0.5, 2, 4)))
    expect_identical(.tile_edges(c(1, 2, 4), visit = TRUE),
        list(lower = c(0.5, 1.5, 3.5), upper = c(1.5, 2.5, 4.5)))
    expect_identical(.tile_edges(5), list(lower = 4.5, upper = 5.5))
})
