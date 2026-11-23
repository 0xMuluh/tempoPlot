# Plot lagged correlations of features with a target

`plotLagHeatmap()` correlates each feature with a target series shifted
by a number of samples within each subject, and draws the correlations
as features by lag tiles.

## Usage

``` r
plotLagHeatmap(x, ...)

# S4 method for class 'SummarizedExperiment'
plotLagHeatmap(
  x,
  target,
  time.col,
  group = NULL,
  assay.type = "counts",
  features = NULL,
  top = 20,
  max.lag = 3,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotLagHeatmap(x, target, ..., experiment = 1L, target.experiment = NULL)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  or
  [`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
  object. For a `MultiAssayExperiment`, `time.col` and `group` are
  columns of its `colData`.

- ...:

  additional arguments passed to
  [`geom_tile`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

- target:

  `Character scalar`. Specifies the target series: a column of
  `colData(x)`, a row of `assay(x, assay.type)`, or
  `"<dimred>:<component>"`.

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects. If `NULL`, all samples form one series. (Default:
  `NULL`)

- assay.type:

  `Character scalar`. Specifies the assay of the features. (Default:
  `"counts"`)

- features:

  `Character vector`. Features to correlate. If `NULL`, uses the `top`
  features with the largest mean. (Default: `NULL`)

- top:

  `Integer scalar`. Number of features when `features` is `NULL`.
  (Default: `20`)

- max.lag:

  `Integer scalar`. Largest lag in both directions, in samples or, for
  visits, in visits. (Default: `3`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

- target.experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment holding `target`, looked up there as for a
  `SummarizedExperiment`, in its first assay for a feature. Only samples
  whose primary sample is in both experiments are used. If `NULL`,
  `target` is looked up in `experiment`. (Default: `NULL`)

## Value

A `ggplot` object.

## Details

At lag \\(k\\), each feature at a sample is paired with the target
\\(k\\) samples later in the same subject. For visits (a factor
`time.col`), it is paired with the target \\(k\\) visits later, and
samples without a sample \\(k\\) visits away are left out. Feature and
target are centred on their mean over the pairs of each subject, so that
differences between subjects do not drive the correlation, and the
Spearman correlation is computed over all pairs. Positive lags show
features that change before the target, negative lags features that
follow it. For other times, lags count samples, so they are comparable
across subjects only for even sampling.

## Examples

``` r
# Species that change before or after acetate in three bioreactors
data("minimalgut", package = "miaTime")
tse <- mia::transformAssay(minimalgut, method = "clr", pseudocount = 1)
plotLagHeatmap(tse, target = "Acetate", time.col = "Time.hr",
    group = "StudyIdentifier", assay.type = "clr")


# The same with acetate from a separate experiment
cd <- SummarizedExperiment::colData(minimalgut)
scfa <- SummarizedExperiment::SummarizedExperiment(
    t(as.matrix(cd[, c("Acetate", "Butyrate", "Propionate")])))
mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(species = tse, SCFA = scfa), colData = cd)
plotLagHeatmap(mae, target = "Acetate", time.col = "Time.hr",
    group = "StudyIdentifier", experiment = "species",
    target.experiment = "SCFA", assay.type = "clr")
```
