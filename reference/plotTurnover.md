# Plot gains and losses of features over time

`plotTurnover()` counts the features gained and lost between consecutive
samples of each subject and draws them as bars spanning the interval
between the samples.

## Usage

``` r
plotTurnover(x, ...)

# S4 method for class 'SummarizedExperiment'
plotTurnover(
  x,
  time.col,
  group = NULL,
  assay.type = "counts",
  detection = 0,
  type = c("gainloss", "components"),
  level = c("subject", "group"),
  group.by = NULL,
  max.gap = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotTurnover(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_rect`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

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

  `Character scalar`. Specifies the assay. (Default: `"counts"`)

- detection:

  `Numeric scalar`. A feature is present in a sample when its value is
  above this. (Default: `0`)

- type:

  `Character scalar`. Draws the numbers of features gained (up) and lost
  (down) (`"gainloss"`), or turnover as the share of features present in
  either sample that appeared and disappeared (`"components"`).
  (Default: `"gainloss"`)

- level:

  `Character scalar`. Draws one panel per subject (`"subject"`) or
  averages the subjects of each `group.by` group at each time
  (`"group"`). (Default: `"subject"`)

- group.by:

  `Character scalar`. Specifies the column from `colData(x)` that groups
  subjects when `level = "group"`. If `NULL`, all subjects form one
  group. (Default: `NULL`)

- max.gap:

  `Numeric scalar`. Skips the intervals between consecutive samples of a
  subject that are further apart than this, in `time.col` units. For
  visits, `NULL` skips intervals across visits without a sample;
  otherwise `NULL` keeps all of them. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

With `type = "components"`, the bars stack appearance and disappearance;
their sum is the total turnover of Hallett et al. (2016).

## References

Hallett LM, et al. (2016). codyn: An R package of community dynamics
metrics. *Methods in Ecology and Evolution* 7(10):1146-1151.

## Examples

``` r
# 16 species in three bioreactors
data("minimalgut", package = "miaTime")
tse <- minimalgut
plotTurnover(tse, time.col = "Time.hr", group = "StudyIdentifier")


# Turnover components, averaged over the bioreactors
plotTurnover(tse, time.col = "Time.hr", group = "StudyIdentifier",
    type = "components", level = "group")
```
