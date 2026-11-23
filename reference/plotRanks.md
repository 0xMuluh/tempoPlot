# Plot the ranks of features over time

`plotRanks()` draws the abundance rank of the most abundant features in
each sample as bump lines, one panel per subject. In polar coordinates
the plot becomes a rank clock.

## Usage

``` r
plotRanks(x, ...)

# S4 method for class 'SummarizedExperiment'
plotRanks(
  x,
  time.col,
  group = NULL,
  assay.type = "counts",
  features = NULL,
  top = 5,
  coord = c("cartesian", "polar"),
  max.gap = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotRanks(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html).

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

  `Character scalar`. Specifies the assay of abundances. To rank taxa
  such as phyla, agglomerate `x` first with
  [`mia::agglomerateByRank()`](https://microbiome.github.io/mia/reference/agglomerate-methods.html).
  (Default: `"counts"`)

- features:

  `Character vector`. Features to draw. If `NULL`, draws the `top`
  features with the best mean rank. (Default: `NULL`)

- top:

  `Integer scalar`. Number of features to draw when `features` is
  `NULL`. (Default: `5`)

- coord:

  `Character scalar`. Draws time on a straight (`"cartesian"`) or
  circular (`"polar"`) axis. (Default: `"cartesian"`)

- max.gap:

  `Numeric scalar`. Breaks a line where consecutive samples of a subject
  are further apart than this, in `time.col` units. For visits, `NULL`
  breaks lines at visits without a sample; otherwise `NULL` never breaks
  lines. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Ranks are among all features of a sample, 1 the most abundant; ties get
their mean rank.

## References

Collins SL, et al. (2008). Rank clocks and plant community dynamics.
*Ecology* 89(12):3534-3541.

## Examples

``` r
# 16 species in three bioreactors
data("minimalgut", package = "miaTime")
tse <- minimalgut
plotRanks(tse, time.col = "Time.hr", group = "StudyIdentifier")


# As rank clocks
plotRanks(tse, time.col = "Time.hr", group = "StudyIdentifier",
    coord = "polar")


# Phyla of the subjects sampled five times
data("hitchip1006", package = "miaTime")
tse <- hitchip1006[, !duplicated(t(SummarizedExperiment::assay(
    hitchip1006, "counts")))]
tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]
tse <- mia::agglomerateByRank(tse, rank = "Phylum")
plotRanks(tse, time.col = "time", group = "subject", top = 3)
```
