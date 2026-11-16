# Extracted from test-plotSwimmer.R:202

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
data("hitchip1006", package = "miaTime")
x <- hitchip1006[, hitchip1006$subject %in%
    names(which(table(hitchip1006$subject) >= 3))]
x$age.time <- x$age + x$time
n_subjects <- length(unique(x$subject))
v <- x[, seq_len(9L)]
v$id <- rep(c("s1", "s2", "s3"), each = 3L)
v$visit <- factor(rep(c("v1", "v2", "v3"), 3L), levels = c("v1", "v2", "v3"))
v <- v[, -c(5L, 9L)]

# test -------------------------------------------------------------------------
p <- plotSwimmer(v, "visit", "id",
        events = data.frame(start = "v2", end = "v3"))
band <- layer_data(p, 1L)
expect_identical(c(band$xmin, band$xmax), c(2L, 3L))
