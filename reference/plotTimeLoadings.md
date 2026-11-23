# Plot temporal loadings of a tensor decomposition

`plotTimeLoadings()` draws the time loading function of each component
of a tensor decomposition stored in `metadata(x)`, such as TEMPTED or
PARAFAC: how strongly the component is expressed at each time.

## Usage

``` r
plotTimeLoadings(x, ...)

# S4 method for class 'SummarizedExperiment'
plotTimeLoadings(x, source, ncomponents = NULL, ...)

# S4 method for class 'MultiAssayExperiment'
plotTimeLoadings(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html).

- source:

  `Character scalar`. Specifies the decomposition in `metadata(x)`: a
  result of
  [`tempted::tempted_all()`](https://rdrr.io/pkg/tempted/man/tempted_all.html)
  or
  [`tempted::tempted()`](https://rdrr.io/pkg/tempted/man/tempted.html),
  a PARAFAC model with the loadings in `Fac` (as from
  `parafac4microbiome::parafac()`), whose third mode is time, or a
  matrix of time points by components with the times as row names or in
  a column `time`.

- ncomponents:

  `Numeric scalar` or `vector`. The components to draw: the first
  `ncomponents`, or the components given by index or name. If `NULL`,
  draws all. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

For TEMPTED, the loadings are functions evaluated on a time grid, and
the legend gives the \\(R^2\\) of each component (`r_square`), of the
data left after the previous components; they do not add up to the total
\\(R^2\\) (`accum_r_square`). For PARAFAC, the time mode has one loading
per time point; set its row names to the times, otherwise points are
numbered in their order. The sign of a component is arbitrary: a
component and its negative describe the same pattern. Feature and
subject loadings are drawn by `miaViz::plotLoadings()` and
[`plotSubjectScatter`](https://0xmuluh.github.io/tempoPlot/reference/plotSubjectScatter.md).

## References

Shi P, et al. (2024). TEMPTED: time-informed dimensionality reduction
for longitudinal microbiome studies. *Genome Biology* 25:317.

## Examples

``` r
# TEMPTED of the species in three bioreactors
if (requireNamespace("tempted", quietly = TRUE)) {
    data("minimalgut", package = "miaTime")
    tse <- minimalgut
    counts <- t(SummarizedExperiment::assay(tse, "counts"))
    S4Vectors::metadata(tse)$tempted <- tempted::tempted_all(counts,
        timepoint = tse$Time_hr_num,
        subjectID = as.character(tse$StudyIdentifier), r = 2,
        threshold = 1, pseudo = 0.5, do_ratio = FALSE)
    plotTimeLoadings(tse, source = "tempted")
}
#> Calculate the 1th Component
#> Convergence reached at dif=8.87596940345795e-08, iter=2
#> Calculate the 2th Component
#> Convergence reached at dif=2.19221530172096e-05, iter=6
```
