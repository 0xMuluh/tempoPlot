library(SummarizedExperiment)
library(ggplot2)

data("minimalgut", package = "miaTime")
x <- minimalgut

test_that("plotReport checks its arguments", {
    skip_if_not_installed("patchwork")
    expect_error(plotReport(x, "Time_hr_num", NULL), "'group'")
    expect_error(plotReport(x, "Time_hr_num", "StudyIdentifier",
        subject = "none"), "'subject'")
    expect_error(plotReport(x, "Time_hr_num", "StudyIdentifier",
        level = "intervention"), "group.by")
})

test_that("plotReport stacks the panels of each level", {
    skip_if_not_installed("patchwork")
    p <- plotReport(x, "Time_hr_num", "StudyIdentifier",
        y = c("Acetate", "Butyrate"), subject = "Bioreactor A")
    expect_s3_class(p, "patchwork")
    expect_length(p$patches$plots, 2L)
    expect_setequal(p$patches$plots[[1L]]$data$StudyIdentifier,
        "Bioreactor A")
    p <- plotReport(x, "Time_hr_num", "StudyIdentifier", y = "Acetate",
        level = "cohort")
    # Swimmer, trend and composition.
    expect_length(p$patches$plots, 2L)
    x$arm <- ifelse(x$StudyIdentifier == "Bioreactor A", "a", "b")
    p <- plotReport(x, "Time_hr_num", "StudyIdentifier", y = "Acetate",
        level = "intervention", group.by = "arm")
    expect_s3_class(p, "patchwork")
})
