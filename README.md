# tempoPlot

<!-- badges: start -->
[![rworkflows](https://github.com/0xMuluh/tempoPlot/actions/workflows/rworkflows.yml/badge.svg?branch=devel)](https://github.com/0xMuluh/tempoPlot/actions)
<!-- badges: end -->

## Microbiome Time Series Visualization

The scope of this package is the visualization of longitudinal microbiome
data. Plots use results computed with
[`miaTime`](https://github.com/microbiome/miaTime) and
[`mia`](https://github.com/microbiome/mia), and the main classes for
interfacing are `TreeSummarizedExperiment` and `MultiAssayExperiment`.

## Installation

```
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("0xMuluh/tempoPlot")
```

## Example

miaTime computes a longitudinal quantity and stores it in the object;
tempoPlot plots it.

```r
library(tempoPlot)
data("minimalgut", package = "miaTime")

# Bray-Curtis dissimilarity of each sample from the previous one
tse <- miaTime::addStepwiseDivergence(minimalgut,
    group = "StudyIdentifier", time.col = "Time.hr",
    assay.type = "counts")

events <- data.frame(start = c(0, 152), end = c(24, NA),
    label = c("batch", "Blautia added"))
plotTrend(tse, y = "divergence", time.col = "Time.hr",
    group = "StudyIdentifier", group.by = "StudyIdentifier",
    events = events)
```

The figure has one line per bioreactor: the dissimilarity of each sample
from the previous one over the 460 h run, between about 0.1 and 1. The
batch phase is shaded, and a vertical line marks the addition of Blautia
hydrogenotrophica at 152 h. The first sample of each bioreactor has no
previous sample and is left out.

The vignette, `vignette("tempoPlot")`, follows a full analysis of a
human cohort and of the bioreactors.

## Which function for which question

| Question | Function |
|---|---|
| When were samples taken? | `plotSwimmer()` |
| | `plotGaps()` |
| | `plotCalendar()` |
| How does a value change over time? | `plotTrend()` |
| | `plotTimeEffect()` |
| How does the composition change? | `plotComposition()` |
| | `plotTiles()` |
| | `plotHorizon()` |
| | `plotRanks()` |
| How fast and how much do communities change? | `plotTurnover()` |
| | `plotTBI()` |
| | `plotLagProfile()` |
| | `plotRecurrence()` |
| | `plotIndividuality()` |
| | `plotFingerprint()` |
| | `plotDOC()` |
| Are there states and rhythms? | `plotTransitionMatrix()` |
| | `plotAlluvial()` |
| | `plotSurvival()` |
| | `plotLandscape()` |
| | `plotSpectrum()` |
| | `plotLagHeatmap()` |
| Which features or subjects stand out? | `plotStats()` |
| | `plotVolcano()` |
| | `plotFeatureScatter()` |
| | `plotSubjectScatter()` |
| | `plotTimeLoadings()` |
| How do subjects move in ordinations and networks? | `plotPath()` |
| | `plotDirection()` |
| | `plotProcrustes()` |
| | `plotNetwork()` |

## Code of conduct

Please note that the tempoPlot project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
