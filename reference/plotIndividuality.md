# Plot within- against between-subject dissimilarities

`plotIndividuality()` compares the dissimilarities between samples of
the same subject with those between samples of different subjects. The
further apart the two distributions, the more individual the
communities.

## Usage

``` r
plotIndividuality(x, ...)

# S4 method for class 'SummarizedExperiment'
plotIndividuality(
  x,
  dissimilarity,
  group,
  layout = c("density", "box", "violin", "time"),
  group.by = NULL,
  time.col = NULL,
  facet.by = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotIndividuality(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_density`](https://ggplot2.tidyverse.org/reference/geom_density.html),
  [`geom_boxplot`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html),
  [`geom_violin`](https://ggplot2.tidyverse.org/reference/geom_violin.html)
  or, for `layout = "time"`,
  [`geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html).

- dissimilarity:

  `Character scalar`. Specifies the sample by sample dissimilarity in
  `metadata(x)`, as stored by
  [`mia::addDissimilarity()`](https://microbiome.github.io/mia/reference/getDissimilarity.html).

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects.

- layout:

  `Character scalar`. Draws density curves (`"density"`), box plots
  (`"box"`), violins (`"violin"`) or the median and interquartile range
  in each time window (`"time"`). (Default: `"density"`)

- group.by:

  `Character scalar`. Specifies the column from `colData(x)` that splits
  the between-subject pairs into pairs of the same and of different
  groups, such as the study centre. (Default: `NULL`)

- time.col:

  `Character scalar`. For `layout = "time"`, specifies the column from
  `colData(x)` holding the time window of each sample: a factor, or
  numeric, `Date` or `POSIXct` values of which each is one window.
  (Default: `NULL`)

- facet.by:

  `Character scalar`. Specifies the column from `colData(x)` that splits
  the pairs into panels. Only pairs of samples from the same group are
  kept. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Each pair of samples is counted once. A subject with \\(k\\) samples
gives \\(k(k-1)/2\\) within-subject pairs, so subjects sampled often
weigh most in the within distribution.

With `layout = "time"`, only pairs of samples from the same window are
kept, and each line is the median dissimilarity of a kind of pair across
the windows, as in Vatanen et al. (2018, doi:10.1038/s41586-018-0620-2,
Fig. 2b). Make windows from times with
[`cut()`](https://rdrr.io/r/base/cut.html).

## Examples

``` r
data("hitchip1006", package = "miaTime")
tse <- hitchip1006
tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 2))]
tse <- mia::addDissimilarity(tse, method = "bray", assay.type = "counts")
plotIndividuality(tse, dissimilarity = "bray", group = "subject")


# By sex
plotIndividuality(tse, dissimilarity = "bray", group = "subject",
    layout = "box", facet.by = "sex")


# Across age windows, with pairs of subjects of the same and of
# different sex
tse$age.window <- cut(tse$age, c(18, 30, 40, 50, 60, 80),
    include.lowest = TRUE)
plotIndividuality(tse, dissimilarity = "bray", group = "subject",
    layout = "time", time.col = "age.window", group.by = "sex")
```
