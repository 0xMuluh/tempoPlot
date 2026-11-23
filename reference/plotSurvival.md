# Plot time to an event

`plotSurvival()` draws Kaplan-Meier curves of the time from each
subject's first sample to an event: dropout, the arrival or loss of a
feature, leaving the first state, or a clinical outcome.

## Usage

``` r
plotSurvival(x, ...)

# S4 method for class 'SummarizedExperiment'
plotSurvival(
  x,
  time.col,
  group,
  event = c("dropout", "arrival", "loss", "state.exit", "outcome"),
  group.by = NULL,
  feature = NULL,
  assay.type = "counts",
  detection = 0,
  state.col = NULL,
  event.col = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotSurvival(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_step`](https://ggplot2.tidyverse.org/reference/geom_path.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects.

- event:

  `Character scalar`. The event: the last sample of a subject is before
  the latest sampling time in `x` (`"dropout"`); `feature` is present
  after being absent at the first sample (`"arrival"`) or absent after
  being present (`"loss"`); `state.col` differs from the first sample
  (`"state.exit"`); or `event.col` is `TRUE` (`"outcome"`). (Default:
  `"dropout"`)

- group.by:

  `Character scalar`. Specifies the column from `colData(x)` that splits
  the subjects into curves, taken from their first sample. (Default:
  `NULL`)

- feature:

  `Character scalar`. Specifies the feature for `event = "arrival"` and
  `"loss"`. A feature is present when its value in `assay.type` is above
  `detection`. (Default: `NULL`)

- assay.type:

  `Character scalar`. Specifies the assay for `feature`. (Default:
  `"counts"`)

- detection:

  `Numeric scalar`. Specifies the value above which `feature` is
  present. (Default: `0`)

- state.col:

  `Character scalar`. Specifies the column from `colData(x)` with the
  state of each sample, for `event = "state.exit"`. (Default: `NULL`)

- event.col:

  `Character scalar`. Specifies the logical column from `colData(x)`
  that marks the samples where the outcome has occurred, for
  `event = "outcome"`. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Follow-up starts at the first sample of each subject. A subject without
the event is censored at its last sample, marked with a cross. The band
is the 95\\ are only seen at sampling times, so the event time is the
first sample at which the event is observed, an upper bound of the true
time.

A subject drops out when its last sample precedes the last sampling time
of the study, the latest time in `x`. A subject that joins late and is
sampled until the end is not a dropout. For `"arrival"` and `"loss"`,
only subjects without and with the feature at their first sample are at
risk. Times are in `time.col` units: days for `Date`, hours for
`POSIXct` times and visit steps for visits.

## References

Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete
observations. *Journal of the American Statistical Association*
53(282):457-481.

## Examples

``` r
# Time until the dominant genus changes, by sex
data("hitchip1006", package = "miaTime")
tse <- mia::addDominant(hitchip1006, name = "dominant")
tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 2))]
plotSurvival(tse, time.col = "time", group = "subject",
    event = "state.exit", state.col = "dominant", group.by = "sex")
```
