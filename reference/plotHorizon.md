# Plot features over time as horizon strips

`plotHorizon()` draws one strip per feature. Each value is centred, cut
into `bands` bands of equal height and the bands are overlaid, so that
deviations above and below the centre show as darker reds and blues in a
short strip.

## Usage

``` r
plotHorizon(x, ...)

# S4 method for class 'SummarizedExperiment'
plotHorizon(
  x,
  time.col,
  group = NULL,
  assay.type = "counts",
  features = NULL,
  top = 20,
  bands = 3,
  centre = c("median", "baseline", "mean"),
  level = c("subject", "group"),
  group.by = NULL,
  max.gap = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotHorizon(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_ribbon`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits. Time differences, as in `max.gap`, are in days for `Date` and
  hours for `POSIXct` times.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects. If `NULL`, all samples form one series. (Default:
  `NULL`)

- assay.type:

  `Character scalar`. Specifies the assay to plot. (Default: `"counts"`)

- features:

  `Character vector`. Features to plot. If `NULL`, plots the `top`
  features with the largest mean. (Default: `NULL`)

- top:

  `Integer scalar`. Number of features to plot when `features` is
  `NULL`. (Default: `20`)

- bands:

  `Integer scalar`. Number of bands on each side of the centre.
  (Default: `3`)

- centre:

  `Character scalar`. Centres each series on its median (`"median"`),
  its first value (`"baseline"`) or its mean (`"mean"`). (Default:
  `"median"`)

- level:

  `Character scalar`. Draws one column of strips per subject
  (`"subject"`) or averages the subjects of each `group.by` group at
  each time (`"group"`). (Default: `"subject"`)

- group.by:

  `Character scalar`. Specifies the column from `colData(x)` that groups
  subjects when `level = "group"`. If `NULL`, all subjects form one
  group. (Default: `NULL`)

- max.gap:

  `Numeric scalar`. Breaks a strip where consecutive samples of a
  subject are further apart than this, in `time.col` units. For visits,
  `NULL` breaks strips at visits without a sample; otherwise `NULL`
  never breaks strips. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

The band height of a feature is its largest absolute deviation from the
centre, over all subjects, divided by `bands`; strips of a feature are
comparable across subjects, not across features. Visits are evenly
spaced.

## References

Fink I, et al. (2022). BiomeHorizon: visualizing microbiome time series
data in R. *mSystems* 7(3):e01380-21.

## Examples

``` r
# 16 species in three bioreactors
data("minimalgut", package = "miaTime")
tse <- mia::transformAssay(minimalgut, method = "relabundance")
plotHorizon(tse, time.col = "Time.hr", group = "StudyIdentifier",
    assay.type = "relabundance", top = 8)


# Deviation from the first sample, averaged over the bioreactors
plotHorizon(tse, time.col = "Time.hr", group = "StudyIdentifier",
    assay.type = "relabundance", top = 8, centre = "baseline",
    level = "group")
```
