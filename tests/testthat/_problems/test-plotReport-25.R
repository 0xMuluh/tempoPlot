# Extracted from test-plotReport.R:25

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
data("minimalgut", package = "miaTime")
x <- minimalgut

# test -------------------------------------------------------------------------
skip_if_not_installed("patchwork")
p <- plotReport(x, "Time_hr_num", "StudyIdentifier",
        y = c("Acetate", "Butyrate"), subject = "Bioreactor A")
expect_s3_class(p, "patchwork")
expect_length(p$patches$plots, 2L)
expect_setequal(p$patches$plots[[1L]]$data$StudyIdentifier,
        "Bioreactor A")
p <- plotReport(x, "Time_hr_num", "StudyIdentifier", y = "Acetate",
        level = "cohort")
