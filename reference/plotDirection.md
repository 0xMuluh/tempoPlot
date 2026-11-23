# Plot directions of movement in an ordination

`plotDirection()` draws a rose diagram of the directions in which
subjects move between consecutive samples in two components of a reduced
dimension.

## Usage

``` r
plotDirection(x, ...)

# S4 method for class 'SummarizedExperiment'
plotDirection(
  x,
  time.col,
  group,
  dimred,
  ncomponents = 2,
  bins = 16,
  facet.by = NULL,
  max.gap = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotDirection(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object, such as a `TreeSummarizedExperiment`.

- ...:

  additional arguments passed to
  [`geom_histogram`](https://ggplot2.tidyverse.org/reference/geom_histogram.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits. Time differences, as in `max.gap`, are in days for `Date` and
  hours for `POSIXct` times.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects.

- dimred:

  `Character scalar`. Specifies the `reducedDim(x)` to use.

- ncomponents:

  `Numeric scalar` or `vector`. The two components of `dimred`: `2` for
  the first two, or two components by index or name. (Default: `2`)

- bins:

  `Integer scalar`. Number of direction sectors. (Default: `16`)

- facet.by:

  `Character scalar`. Specifies the column from `colData(x)` that splits
  the moves into panels. Each move takes the value of its later sample.
  (Default: `NULL`)

- max.gap:

  `Numeric scalar`. Consecutive samples of a subject further apart than
  this, in `time.col` units, are not a move. For visits, `NULL` skips
  moves across visits without a sample; otherwise `NULL` counts all.
  (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Angles are measured counterclockwise from the first component (east), so
90 degrees points along the second component. Moves of zero length have
no direction and are left out. The radius is on a square-root scale, so
that the area of a sector follows its number of moves.

## Examples

``` r
# Directions in a PCoA of three bioreactors
data("minimalgut", package = "miaTime")
tse <- mia::transformAssay(minimalgut, method = "relabundance")
tse <- mia::addMDS(tse, assay.type = "relabundance", method = "bray",
    name = "PCoA")
plotDirection(tse, time.col = "Time.hr", group = "StudyIdentifier",
    dimred = "PCoA", facet.by = "StudyIdentifier")
```
