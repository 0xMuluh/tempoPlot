library(SummarizedExperiment)
library(ggplot2)

# Subjects f1 and f2, sampled at times 1 and 2, with partners p1 (one
# sample) and p2 (two samples), whose samples have no time.
ids <- c("f1a", "f1b", "f2a", "f2b", "p1", "p2a", "p2b")
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 7L)),
    colData = DataFrame(id = c("f1", "f1", "f2", "f2", "p1", "p2", "p2"),
        t = c(1, 2, 1, 2, NA, NA, NA),
        partner = c("p1", "p1", "p2", "p2", NA, NA, NA),
        site = c("x", "x", "y", "y", NA, NA, NA)))
colnames(toy) <- ids
d <- matrix(0.5, 7L, 7L, dimnames = list(ids, ids))
diag(d) <- 0
# Own partner: means 0.2, 0.3, 0.2, 0.4. Other partners: 0.7, 0.7, 0.9,
# 0.5.
pairs <- rbind(c("f1a", "p1", 0.2), c("f1a", "p2a", 0.6),
    c("f1a", "p2b", 0.8), c("f1b", "p1", 0.3), c("f1b", "p2a", 0.5),
    c("f1b", "p2b", 0.9), c("f2a", "p2a", 0.1), c("f2a", "p2b", 0.3),
    c("f2a", "p1", 0.9), c("f2b", "p2a", 0.4), c("f2b", "p2b", 0.4),
    c("f2b", "p1", 0.5))
d[pairs[, 1:2]] <- d[pairs[, 2:1]] <- as.numeric(pairs[, 3])
S4Vectors::metadata(toy)$d <- d

test_that("plotPartner checks its arguments", {
    expect_error(plotPartner(toy, "none", "t", "id", "partner"),
        "'dissimilarity'")
    expect_error(plotPartner(toy, "d", "t", "id", "none"), "'partner.col'")
    expect_error(plotPartner(toy, "d", "t", "id", "partner",
        layout = "bar"))
    bad <- toy
    bad$partner[2L] <- "p2"
    expect_error(plotPartner(bad, "d", "t", "id", "partner"),
        "constant within each subject")
    bad$partner[1:2] <- "p3"
    expect_error(plotPartner(bad, "d", "t", "id", "partner"),
        "subjects of the 'group' column")
    bad$partner <- NA
    expect_error(plotPartner(bad, "d", "t", "id", "partner"),
        "at least one subject")
})

test_that("plotPartner averages over own and other partners' samples", {
    # Partners have a single sample or no time, without a warning.
    expect_no_warning(p <- plotPartner(toy, "d", "t", "id", "partner"))
    df <- p$data
    expect_identical(levels(df$.partner),
        c("own partner", "other partners"))
    own <- df[df$.partner == "own partner", ]
    other <- df[df$.partner == "other partners", ]
    expect_identical(own$id, c("f1", "f1", "f2", "f2"))
    expect_equal(own$t, c(1, 2, 1, 2))
    expect_equal(own$.y, c(0.2, 0.3, 0.2, 0.4))
    expect_equal(other$.y, c(0.7, 0.7, 0.9, 0.5))
    # One line per subject and comparison, in two colours.
    ld <- layer_data(p, 1L)
    expect_length(unique(ld$group), 4L)
    expect_setequal(ld$colour, c("#0072B2", "#E69F00"))
    expect_identical(p$labels$colour, "partner")
    expect_identical(p$labels$y, "Mean d")
})

test_that("plotPartner draws boxes per time and comparison", {
    p <- plotPartner(toy, "d", "t", "id", "partner", layout = "box",
        facet.by = "site")
    expect_s3_class(p$layers[[1L]]$geom, "GeomBoxplot")
    ld <- layer_data(p, 1L)
    expect_identical(nrow(ld), 8L)
    expect_setequal(as.character(p$data$site), c("x", "y"))
})

test_that("plotPartner draws own partners alone without other partners", {
    one <- toy
    one$partner[3:4] <- "p1"
    df <- plotPartner(one, "d", "t", "id", "partner")$data
    expect_identical(as.character(unique(df$.partner)), "own partner")
    expect_equal(df$.y, c(0.2, 0.3, 0.9, 0.5))
})
