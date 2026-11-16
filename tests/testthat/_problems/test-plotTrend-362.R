# Extracted from test-plotTrend.R:362

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)
data("minimalgut", package = "miaTime")
x <- minimalgut
n_subjects <- length(unique(x$StudyIdentifier))
n_times <- length(unique(x$Time_hr_num))

# test -------------------------------------------------------------------------
data("minimalgut", package = "miaTime")
ref <- data.frame(Time_hr_num = c(0, 460), estimate = c(10, 20),
        lower = c(5, 15), upper = c(15, 25))
p <- plotTrend(minimalgut, "Acetate", "Time_hr_num", "StudyIdentifier",
        reference = ref, limits = c(30, 60))
geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
expect_identical(geoms[1:4],
        c("GeomRect", "GeomRibbon", "GeomLine", "GeomLine"))
