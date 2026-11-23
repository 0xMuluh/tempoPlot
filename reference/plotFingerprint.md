# Plot whether subjects are identified by their own earlier sample

`plotFingerprint()` compares, for each subject, the sample at a query
time with all samples at a reference time. A subject is identified when
its own reference sample is closer than the reference sample of any
other subject.

## Usage

``` r
plotFingerprint(x, ...)

# S4 method for class 'SummarizedExperiment'
plotFingerprint(x, dissimilarity, time.col, group, time.point, ...)

# S4 method for class 'MultiAssayExperiment'
plotFingerprint(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

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

- time.point:

  `Vector` of two values of `time.col`: the reference time, whose
  samples subjects are identified from, and the query time, whose
  samples are identified.

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Each point is a subject with samples at both times: the dissimilarity of
its query sample to its own reference sample (x) and to the nearest
reference sample of another subject (y). Points above the 1:1 line are
identified. The subtitle gives the share of subjects identified.
Franzosa et al. (2015) identified people by metagenomic codes; this plot
applies the same idea to the nearest sample on a dissimilarity.

## References

Franzosa EA, et al. (2015). Identifying personal microbiomes using
metagenomic codes. *PNAS* 112(22):E2930-E2938.

## Examples

``` r
# Is each subject's second sample closest to its own first sample?
data("hitchip1006", package = "miaTime")
tse <- hitchip1006
# Some profiles are stored twice; keep the first copy
tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
tse <- tse[, tse$visit %in% c("1", "2")]
tse <- mia::addDissimilarity(tse, method = "bray", assay.type = "counts")
plotFingerprint(tse, dissimilarity = "bray", time.col = "visit",
    group = "subject", time.point = c("1", "2"))
#> Warning: 'group' has 929 subjects with a single sample.
```
