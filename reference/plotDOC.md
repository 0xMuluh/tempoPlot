# Plot dissimilarity against overlap of sample pairs

`plotDOC()` draws the dissimilarity of each pair of samples against
their overlap, with a smooth: the dissimilarity-overlap curve (DOC). A
falling curve at high overlap suggests that the same taxa behave alike
in different communities.

## Usage

``` r
plotDOC(x, ...)

# S4 method for class 'SummarizedExperiment'
plotDOC(
  x,
  overlap,
  dissimilarity,
  group = NULL,
  pairs = c("all", "within", "between"),
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotDOC(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

- overlap, dissimilarity:

  `Character scalar`. Specify the sample by sample matrices in
  `metadata(x)`, as stored by
  [`mia::addDissimilarity()`](https://microbiome.github.io/mia/reference/getDissimilarity.html):
  the overlap on the x axis and the dissimilarity on the y axis.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects, needed for `pairs = "within"` and `"between"`.
  (Default: `NULL`)

- pairs:

  `Character scalar`. Draws all pairs of samples (`"all"`), pairs from
  the same subject (`"within"`) or from different subjects
  (`"between"`). (Default: `"all"`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Each pair of samples is drawn once, with the robust LOWESS curve of
Bashan et al. (2016); above 1000 pairs, points are drawn with
transparency. In Bashan et al. (2016), overlap is mia's `"overlap"` of
relative abundances and dissimilarity is the root Jensen-Shannon
divergence of the shared taxa, renormalised; mia's `"jsd"` uses all
taxa.

## References

Bashan A, et al. (2016). Universality of human microbial dynamics.
*Nature* 534(7606):259-262.

## Examples

``` r
# Samples of four artificial guts
data("SilvermanAGutData", package = "miaTime")
tse <- SilvermanAGutData
# Five samples have no reads
tse <- tse[, colSums(SummarizedExperiment::assay(tse)) > 0]
tse <- mia::transformAssay(tse, method = "relabundance")
tse <- mia::addDissimilarity(tse, method = "overlap",
    assay.type = "relabundance")
tse <- mia::addDissimilarity(tse, method = "jsd",
    assay.type = "relabundance")
plotDOC(tse, overlap = "overlap", dissimilarity = "jsd",
    group = "Vessel", pairs = "between")
```
