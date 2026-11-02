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

# Code of conduct

Please note that the tempoPlot project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
