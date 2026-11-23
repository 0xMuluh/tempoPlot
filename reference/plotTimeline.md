# Plot the timeline of one subject

`plotTimeline()` stacks tracks of one subject on a shared time axis: the
composition of the most abundant features, numeric variables as lines, a
strip of states and a strip of events.

## Usage

``` r
plotTimeline(x, ...)

# S4 method for class 'SummarizedExperiment'
plotTimeline(
  x,
  time.col,
  group = NULL,
  subject = NULL,
  assay.type = "counts",
  y = NULL,
  state.col = NULL,
  events = NULL,
  max.gap = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotTimeline(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`plotComposition`](https://0xmuluh.github.io/tempoPlot/reference/plotComposition.md),
  such as `features`, `top`, `as.relative` or `layout`.

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits. Time differences, as in `max.gap`, are in days for `Date` and
  hours for `POSIXct` times.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects. If `NULL`, all samples belong to one subject.
  (Default: `NULL`)

- subject:

  `Character scalar`. The subject to draw, a value of the `group`
  column. (Default: `NULL`)

- assay.type:

  `Character scalar`. Specifies the assay of abundances, also read for
  `y` that are features. (Default: `"counts"`)

- y:

  `Character vector`. Specifies the variables drawn as lines, one track
  each, as in
  [`plotTrend`](https://0xmuluh.github.io/tempoPlot/reference/plotTrend.md):
  columns of `colData(x)`, rows of `assay(x, assay.type)` or
  `"<dimred>:<component>"`. (Default: `NULL`)

- state.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains the state of each sample, such as a community type, drawn as
  a strip of tiles. (Default: `NULL`)

- events:

  `data.frame`. Events to draw in a strip at the bottom, one per row,
  with columns `start`, and optionally `end`, `label` and the `group`
  column. `start` and `end` have the type of `time.col`; for visits,
  they are visit levels. Events with an `end` are drawn as bands, the
  others as vertical lines; `label` sets their fill or line type. Events
  with a `group` column are drawn for their subject only. (Default:
  `NULL`)

- max.gap:

  `Numeric scalar`. Breaks lines, and areas of the composition, where
  consecutive samples are further apart than this, in `time.col` units.
  For visits, `NULL` breaks them at visits without a sample; otherwise
  `NULL` never breaks them. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `patchwork` object, not a single `ggplot`: the tracks from top to
bottom.

## Details

The tracks share the x axis: the same limits, from half the spacing of
samples before the first to half after the last, widened to include the
events, and the same breaks. Only the bottom track shows the time axis,
and the legends of all tracks are collected on the right. The layout
follows the subject timelines of Gajer et al. (2012), Figs 3 and 5, Zhou
et al. (2019), Fig. 6b, and David et al. (2014), Fig. 1.

Draws with the patchwork package. A single track is changed with
`p[[i]] <- p[[i]] + ...`, all of them with `p & ...`.

## References

David LA, et al. (2014). Host lifestyle affects human microbiota on
daily timescales. *Genome Biology* 15(7):R89.
doi:10.1186/gb-2014-15-7-r89

Gajer P, et al. (2012). Temporal dynamics of the human vaginal
microbiota. *Science Translational Medicine* 4(132):132ra52.
doi:10.1126/scitranslmed.3003605

Zhou W, et al. (2019). Longitudinal multi-omics of host-microbe dynamics
in prediabetes. *Nature* 569(7758):663-671.
doi:10.1038/s41586-019-1236-x

## Examples

``` r
# Daily samples of one woman over six weeks: genera, cell count per
# gram of stool and enterotype
data("temporalMicrobiome20", package = "miaTime")
tse <- temporalMicrobiome20
# Samples without a profile have the enterotype "NA"
tse <- tse[, tse$enterotype != "NA"]
plotTimeline(tse, time.col = "day", group = "person_id", subject = 820,
    y = "cell_count_per_gram", state.col = "enterotype", top = 6)


# A bioreactor of a defined community: species, butyrate, the dominant
# species of each sample and the perturbations of its schedule
data("minimalgut", package = "miaTime")
tse <- mia::addDominant(minimalgut, name = "dominant")
events <- data.frame(
    start = c(0, 152, 432),
    end = c(24, NA, 436),
    label = c("batch", "Blautia added", "double dilution rate"))
plotTimeline(tse, time.col = "Time.hr", group = "StudyIdentifier",
    subject = "Bioreactor A", y = "Butyrate", state.col = "dominant",
    events = events, top = 8, layout = "area")
```
