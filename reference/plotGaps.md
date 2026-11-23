# Plot intervals between consecutive samples

`plotGaps()` draws the distribution of the intervals between consecutive
samples of each subject.

## Usage

``` r
plotGaps(x, ...)

# S4 method for class 'SummarizedExperiment'
plotGaps(
  x,
  time.col,
  group,
  group.by = NULL,
  layout = c("histogram", "density", "box"),
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotGaps(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_histogram`](https://ggplot2.tidyverse.org/reference/geom_histogram.html),
  [`geom_density`](https://ggplot2.tidyverse.org/reference/geom_density.html)
  or
  [`geom_boxplot`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects.

- group.by:

  `Character scalar`. Specifies the column from `colData(x)` that splits
  the intervals into groups. Each interval takes the value of its later
  sample. (Default: `NULL`)

- layout:

  `Character scalar`. Draws a histogram (`"histogram"`), density curves
  (`"density"`) or box plots (`"box"`). (Default: `"histogram"`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Intervals are in `time.col` units: days for `Date`, hours for `POSIXct`
times and visit steps for visits.

## Examples

``` r
data("hitchip1006", package = "miaTime")
tse <- hitchip1006
tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 2))]
plotGaps(tse, time.col = "time", group = "subject", bins = 30)


# Intervals by sex
plotGaps(tse, time.col = "time", group = "subject", group.by = "sex",
    layout = "box")
```
