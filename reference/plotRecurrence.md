# Plot dissimilarities between the samples of each subject

`plotRecurrence()` draws, for each subject, a time by time matrix of the
dissimilarities between its samples. With `threshold`, it draws a
recurrence plot: the pairs of samples closer than the threshold.

## Usage

``` r
plotRecurrence(x, ...)

# S4 method for class 'SummarizedExperiment'
plotRecurrence(
  x,
  dissimilarity,
  time.col,
  group,
  subjects = NULL,
  threshold = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotRecurrence(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_rect`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

- dissimilarity:

  `Character scalar`. Specifies the sample by sample dissimilarity in
  `metadata(x)`, as stored by
  [`mia::addDissimilarity()`](https://microbiome.github.io/mia/reference/getDissimilarity.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects.

- subjects:

  `Character vector`. Subjects to draw, one panel each. If `NULL`, draws
  all. (Default: `NULL`)

- threshold:

  `Numeric scalar`. Marks the pairs of samples with a dissimilarity of
  at most this. If `NULL`, colours by dissimilarity. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Each sample covers the time from halfway to its previous sample to
halfway to its next one, so the axes keep the spacing of the sampling
times. Visits are evenly spaced.

## Examples

``` r
# Three bioreactors under the same perturbations
data("minimalgut", package = "miaTime")
tse <- mia::addDissimilarity(minimalgut, method = "bray",
    assay.type = "counts")
plotRecurrence(tse, dissimilarity = "bray", time.col = "Time.hr",
    group = "StudyIdentifier")


# Pairs of samples within a Bray-Curtis dissimilarity of 0.2
plotRecurrence(tse, dissimilarity = "bray", time.col = "Time.hr",
    group = "StudyIdentifier", subjects = "Bioreactor A",
    threshold = 0.2)
```
