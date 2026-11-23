# Plot flows of subjects between states over time

`plotAlluvial()` draws the state of each subject at each time point as
stacked strata and the subjects moving between states as flows.

## Usage

``` r
plotAlluvial(x, ...)

# S4 method for class 'SummarizedExperiment'
plotAlluvial(
  x,
  state.col,
  time.col,
  group,
  min.flow = 1,
  missing = c("show", "drop"),
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotAlluvial(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`ggalluvial::geom_flow`](http://corybrunson.github.io/ggalluvial/reference/geom_flow.md).

- state.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains the state of each sample, such as a community type or the
  dominant taxon.

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains the time points shared by the subjects: numeric, `Date`,
  `POSIXct` or, best, a factor of visits. Each distinct time is one
  axis.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects.

- min.flow:

  `Integer scalar`. Hides the paths through all time points followed by
  fewer subjects than this. (Default: `1`)

- missing:

  `Character scalar`. Draws time points without a sample or a state as a
  `"missing"` stratum (`"show"`) or removes the subjects that have them
  (`"drop"`). (Default: `"show"`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Requires the ggalluvial package. Subjects that follow the same states
through all time points form one flow.

## Examples

``` r
# Dominant genus at the first three samples of each subject
data("hitchip1006", package = "miaTime")
tse <- mia::addDominant(hitchip1006, name = "dominant", n = 4)
tse$visit <- ave(tse$time, tse$subject, FUN = rank)
tse <- tse[, tse$visit <= 3 &
    tse$subject %in% names(which(table(tse$subject) >= 2))]
tse$visit <- factor(tse$visit)
plotAlluvial(tse, state.col = "dominant", time.col = "visit",
    group = "subject")
```
