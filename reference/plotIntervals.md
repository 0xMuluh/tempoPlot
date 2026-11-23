# Plot time intervals of features

`plotIntervals()` draws time intervals of each feature as horizontal
bars, such as the intervals in which a feature differs between groups.

## Usage

``` r
plotIntervals(x, ...)

# S4 method for class 'SummarizedExperiment'
plotIntervals(x, name, ...)

# S4 method for class 'data.frame'
plotIntervals(
  x,
  feature.col = "feature",
  colour.by = NULL,
  facet.by = NULL,
  features = NULL,
  top = 20,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotIntervals(x, ..., experiment = 1L)
```

## Arguments

- x:

  a `data.frame` of intervals, one row per interval with columns `start`
  and `end`, or a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object with such a table in its `metadata`.

- ...:

  additional arguments passed to
  [`geom_segment`](https://ggplot2.tidyverse.org/reference/geom_segment.html).

- name:

  `Character scalar`. Specifies the entry of `metadata(x)` holding the
  intervals.

- feature.col:

  `Character scalar`. Specifies the column of feature names. (Default:
  `"feature"`)

- colour.by:

  `Character scalar`. Specifies the column that colours the intervals,
  such as the group with the higher values. (Default: `NULL`)

- facet.by:

  `Character scalar`. Specifies the column that splits the intervals
  into panels, such as the method that found them. (Default: `NULL`)

- features:

  `Character vector`. Features to plot. If `NULL`, plots the `top`
  features with the longest total time in intervals. (Default: `NULL`)

- top:

  `Integer scalar`. Number of features to plot when `features` is
  `NULL`. (Default: `20`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

`start` and `end` are numeric, `Date`, `POSIXct` or factors of visits
with the same levels. An interval with `end` equal to `start`, or
missing, covers one time point and is drawn as a point. Features are
ordered by the start of their first interval, the earliest at the top.

The columns `feature`, `start`, `end` and `dominant` match the
significant intervals returned by MetaLonDA (Metwally et al. 2018,
doi:10.1186/s40168-018-0402-y); intervals from SplinectomeR or from
comparisons of fitted curves can be given in the same form.

## Examples

``` r
# Illustrative significant intervals in the form returned by MetaLonDA:
# days in which each genus differs between two groups, and the group
# with the higher abundance
intervals <- data.frame(
    feature = c("Bacteroides", "Bacteroides", "Prevotella",
        "Akkermansia", "Bifidobacterium"),
    start = c(10, 60, 30, 0, 45),
    end = c(25, 90, 75, 20, 45),
    dominant = c("case", "control", "case", "control", "case"))
plotIntervals(intervals, colour.by = "dominant")


# The same intervals stored in metadata(tse)
data("minimalgut", package = "miaTime")
tse <- minimalgut
S4Vectors::metadata(tse)$intervals <- intervals
plotIntervals(tse, name = "intervals", colour.by = "dominant")
```
