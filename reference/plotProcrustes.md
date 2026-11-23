# Plot two ordinations superimposed by Procrustes rotation

`plotProcrustes()` rotates a second ordination of the same samples onto
the first and draws both, with a segment from each sample in the first
to the same sample in the second. Short segments mean the two
ordinations agree.

## Usage

``` r
plotProcrustes(x, ...)

# S4 method for class 'SummarizedExperiment'
plotProcrustes(x, dimred, dimred2, ncomponents = 2, colour.by = NULL, ...)

# S4 method for class 'MultiAssayExperiment'
plotProcrustes(
  x,
  experiment,
  experiment2,
  dimred,
  dimred2 = dimred,
  ncomponents = 2,
  colour.by = NULL,
  ...
)
```

## Arguments

- x:

  a
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object, such as a `TreeSummarizedExperiment`, or a
  [`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html).

- ...:

  additional arguments passed to
  [`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

- dimred:

  `Character scalar`. Specifies the first `reducedDim`.

- dimred2:

  `Character scalar`. Specifies the second `reducedDim`, of `x` or of
  `experiment2`. For a `MultiAssayExperiment`, it defaults to `dimred`.

- ncomponents:

  `Numeric scalar` or `vector`. The components of each ordination used
  in the rotation: the first `ncomponents`, or the components given by
  index or name. The first two are drawn. (Default: `2`)

- colour.by:

  `Character scalar`. Specifies the column from `colData(x)` that
  colours the points. (Default: `NULL`)

- experiment, experiment2:

  `Character scalar`. For a `MultiAssayExperiment`, the experiments
  holding `dimred` and `dimred2`. Samples are matched by their primary
  names.

## Value

A `ggplot` object.

## Details

Both ordinations are centred and scaled to unit sum of squares, and the
second is rotated to the first (symmetric Procrustes analysis). The
subtitle gives the Procrustes statistic \\(m^2\\), the sum of squared
distances between matched samples, from 0 (same configuration) to 1.

## References

Peres-Neto PR, Jackson DA (2001). How well do multivariate data sets
match? The advantages of a Procrustean superimposition approach over the
Mantel test. *Oecologia* 129(2):169-178.

## Examples

``` r
# PCoA of Bray-Curtis against PCA of clr abundances
data("minimalgut", package = "miaTime")
tse <- mia::transformAssay(minimalgut, method = "relabundance")
tse <- mia::transformAssay(tse, method = "clr", pseudocount = 1,
    assay.type = "counts")
tse <- mia::addMDS(tse, assay.type = "relabundance", method = "bray",
    name = "PCoA")
tse <- scater::runPCA(tse, assay.type = "clr", ncomponents = 5,
    name = "PCA")
plotProcrustes(tse, dimred = "PCoA", dimred2 = "PCA",
    colour.by = "StudyIdentifier")
```
