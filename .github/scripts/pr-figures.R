# Run the examples of the given help pages and save each plot as a PNG.
# Usage: Rscript pr-figures.R <output dir> <topic> [<topic> ...]
args <- commandArgs(trailingOnly = TRUE)
out <- args[1L]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
library(tempoPlot)

for (topic in args[-1L]) {
    ex <- tempfile(fileext = ".R")
    tools::Rd2ex(file.path("man", paste0(topic, ".Rd")), ex)
    if (!file.exists(ex)) next
    png(file.path(out, paste0(topic, "-%02d.png")), width = 1600,
        height = 1000, res = 200)
    source(ex, print.eval = TRUE, local = new.env())
    dev.off()
}
