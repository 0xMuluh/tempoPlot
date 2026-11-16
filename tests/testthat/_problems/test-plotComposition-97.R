# Extracted from test-plotComposition.R:97

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
data("minimalgut", package = "miaTime")
mg <- minimalgut
prop <- prop.table(assay(mg, "counts"), 2L)
data("hitchip1006", package = "miaTime")
h <- hitchip1006[, !duplicated(t(assay(hitchip1006, "counts")))]
h <- h[, h$subject %in% names(which(table(h$subject) == 5))]
h$visit <- factor(ave(h$time, h$subject, FUN = rank))

# test -------------------------------------------------------------------------
p <- plotComposition(x, "Time_hr_num", "StudyIdentifier", top = 3,
        layout = "stream")
