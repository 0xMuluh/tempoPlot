# Plot effects against their significance

`plotVolcano()` draws the effect of each feature against
\\(-\log\_{10}\\) of its p-value and labels the strongest hits.

## Usage

``` r
plotVolcano(x, ...)

# S4 method for class 'SummarizedExperiment'
plotVolcano(
  x,
  effect,
  p.col,
  threshold = 0.05,
  label.top = 10,
  feature.col = NULL,
  ...
)

# S4 method for class 'data.frame'
plotVolcano(
  x,
  effect,
  p.col,
  threshold = 0.05,
  label.top = 10,
  feature.col = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotVolcano(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object, with the results in `rowData(x)`, or a `data.frame` of results
  with one row per feature.

- ...:

  additional arguments passed to
  [`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

- effect:

  `Character scalar`. Specifies the column holding the effects, such as
  log fold changes or slopes.

- p.col:

  `Character scalar`. Specifies the column holding the p-values or
  adjusted p-values.

- threshold:

  `Numeric scalar`. Features with `p.col` below this are significant.
  (Default: `0.05`)

- label.top:

  `Integer scalar`. Labels this many significant features with the
  smallest p-values. (Default: `10`)

- feature.col:

  `Character scalar`. Specifies the column of feature names that label
  the points, of `rowData(x)` or of the `data.frame`. If `NULL`, uses
  the row names. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Features with a missing effect or p-value are dropped. The dashed line
marks `threshold`.

## Examples

``` r
# Differences between men and women at the first sample of each subject
data("hitchip1006", package = "miaTime")
tse <- hitchip1006[, hitchip1006$time == 0 & !is.na(hitchip1006$sex)]
tse <- mia::transformAssay(tse, method = "clr", pseudocount = 1)
mat <- SummarizedExperiment::assay(tse, "clr")
male <- tse$sex == "male"
res <- data.frame(
    genus = rownames(tse),
    difference = rowMeans(mat[, male]) - rowMeans(mat[, !male]),
    p = apply(mat, 1, function(v) wilcox.test(v[male], v[!male])$p.value))
res$q <- p.adjust(res$p, method = "fdr")
plotVolcano(res, effect = "difference", p.col = "q", label.top = 5,
    feature.col = "genus")
```
