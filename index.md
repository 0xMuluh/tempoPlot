# tempoPlot

## Microbiome Time Series Visualization

The scope of this package is the visualization of longitudinal
microbiome data. Plots use results computed with
[`miaTime`](https://github.com/microbiome/miaTime) and
[`mia`](https://github.com/microbiome/mia), and the main classes for
interfacing are `TreeSummarizedExperiment` and `MultiAssayExperiment`.

## Installation

    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")

    remotes::install_github("0xMuluh/tempoPlot")

## Example

miaTime computes a longitudinal quantity and stores it in the object;
tempoPlot plots it.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tempoPlot`](https://github.com/0xMuluh/tempoPlot)`)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"minimalgut"``, package ``=`` ``"miaTime"``)`\
\
`# Bray-Curtis dissimilarity of each sample from the previous one`\
`tse`` ``<-`` ``miaTime``::`[`addStepwiseDivergence`](https://rdrr.io/pkg/miaTime/man/getStepwiseDivergence.html)`(``minimalgut``,`\
`    group ``=`` ``"StudyIdentifier"``, time.col ``=`` ``"Time.hr"``,`\
`    assay.type ``=`` ``"counts"``)`\
\
`events`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``start ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``152``)``, end ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``24``, ``NA``)``,`\
`    label ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"batch"``, ``"Blautia added"``)``)`\
[`plotTrend`](https://0xmuluh.github.io/tempoPlot/reference/plotTrend.md)`(``tse``, y ``=`` ``"divergence"``, time.col ``=`` ``"Time.hr"``,`\
`    group ``=`` ``"StudyIdentifier"``, group.by ``=`` ``"StudyIdentifier"``,`\
`    events ``=`` ``events``)`

The figure has one line per bioreactor: the dissimilarity of each sample
from the previous one over the 460 h run, between about 0.1 and 1. The
batch phase is shaded, and a vertical line marks the addition of Blautia
hydrogenotrophica at 152 h. The first sample of each bioreactor has no
previous sample and is left out.

The vignette,
[`vignette("tempoPlot")`](https://0xmuluh.github.io/tempoPlot/articles/tempoPlot.md),
follows a full analysis of a human cohort and of the bioreactors.

## Which function for which question

| Question | Function |
|----|----|
| When were samples taken? | [`plotSwimmer()`](https://0xmuluh.github.io/tempoPlot/reference/plotSwimmer.md) |
|  | [`plotGaps()`](https://0xmuluh.github.io/tempoPlot/reference/plotGaps.md) |
|  | [`plotCalendar()`](https://0xmuluh.github.io/tempoPlot/reference/plotCalendar.md) |
| What happened to one subject? | [`plotTimeline()`](https://0xmuluh.github.io/tempoPlot/reference/plotTimeline.md) |
| How does a value change over time? | [`plotTrend()`](https://0xmuluh.github.io/tempoPlot/reference/plotTrend.md) |
|  | [`plotTimeEffect()`](https://0xmuluh.github.io/tempoPlot/reference/plotTimeEffect.md) |
|  | [`plotIntervals()`](https://0xmuluh.github.io/tempoPlot/reference/plotIntervals.md) |
|  | [`plotRmcorr()`](https://0xmuluh.github.io/tempoPlot/reference/plotRmcorr.md) |
| How does the composition change? | [`plotComposition()`](https://0xmuluh.github.io/tempoPlot/reference/plotComposition.md) |
|  | [`plotTiles()`](https://0xmuluh.github.io/tempoPlot/reference/plotTiles.md) |
|  | [`plotProfiles()`](https://0xmuluh.github.io/tempoPlot/reference/plotProfiles.md) |
|  | [`plotHorizon()`](https://0xmuluh.github.io/tempoPlot/reference/plotHorizon.md) |
|  | [`plotRanks()`](https://0xmuluh.github.io/tempoPlot/reference/plotRanks.md) |
| How fast and how much do communities change? | [`plotTurnover()`](https://0xmuluh.github.io/tempoPlot/reference/plotTurnover.md) |
|  | [`plotChanges()`](https://0xmuluh.github.io/tempoPlot/reference/plotChanges.md) |
|  | [`plotTBI()`](https://0xmuluh.github.io/tempoPlot/reference/plotTBI.md) |
|  | [`plotLagProfile()`](https://0xmuluh.github.io/tempoPlot/reference/plotLagProfile.md) |
|  | [`plotRecurrence()`](https://0xmuluh.github.io/tempoPlot/reference/plotRecurrence.md) |
|  | [`plotIndividuality()`](https://0xmuluh.github.io/tempoPlot/reference/plotIndividuality.md) |
|  | [`plotFingerprint()`](https://0xmuluh.github.io/tempoPlot/reference/plotFingerprint.md) |
|  | [`plotDOC()`](https://0xmuluh.github.io/tempoPlot/reference/plotDOC.md) |
| Are there states and rhythms? | [`plotTransitionMatrix()`](https://0xmuluh.github.io/tempoPlot/reference/plotTransitionMatrix.md) |
|  | [`plotAlluvial()`](https://0xmuluh.github.io/tempoPlot/reference/plotAlluvial.md) |
|  | [`plotStates()`](https://0xmuluh.github.io/tempoPlot/reference/plotStates.md) |
|  | [`plotRuns()`](https://0xmuluh.github.io/tempoPlot/reference/plotRuns.md) |
|  | [`plotSurvival()`](https://0xmuluh.github.io/tempoPlot/reference/plotSurvival.md) |
|  | [`plotLandscape()`](https://0xmuluh.github.io/tempoPlot/reference/plotLandscape.md) |
|  | [`plotSpectrum()`](https://0xmuluh.github.io/tempoPlot/reference/plotSpectrum.md) |
|  | [`plotNoise()`](https://0xmuluh.github.io/tempoPlot/reference/plotNoise.md) |
|  | [`plotLagHeatmap()`](https://0xmuluh.github.io/tempoPlot/reference/plotLagHeatmap.md) |
| What does a subject share with its partner? | [`plotPartner()`](https://0xmuluh.github.io/tempoPlot/reference/plotPartner.md) |
|  | [`plotOrigin()`](https://0xmuluh.github.io/tempoPlot/reference/plotOrigin.md) |
| Which features or subjects stand out? | [`plotStats()`](https://0xmuluh.github.io/tempoPlot/reference/plotStats.md) |
|  | [`plotVolcano()`](https://0xmuluh.github.io/tempoPlot/reference/plotVolcano.md) |
|  | [`plotFeatureScatter()`](https://0xmuluh.github.io/tempoPlot/reference/plotFeatureScatter.md) |
|  | [`plotSubjectScatter()`](https://0xmuluh.github.io/tempoPlot/reference/plotSubjectScatter.md) |
|  | [`plotTimeLoadings()`](https://0xmuluh.github.io/tempoPlot/reference/plotTimeLoadings.md) |
| How do subjects move in ordinations and networks? | [`plotPath()`](https://0xmuluh.github.io/tempoPlot/reference/plotPath.md) |
|  | [`plotDirection()`](https://0xmuluh.github.io/tempoPlot/reference/plotDirection.md) |
|  | [`plotProcrustes()`](https://0xmuluh.github.io/tempoPlot/reference/plotProcrustes.md) |
|  | [`plotNetwork()`](https://0xmuluh.github.io/tempoPlot/reference/plotNetwork.md) |

## Code of conduct

Please note that the tempoPlot project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
