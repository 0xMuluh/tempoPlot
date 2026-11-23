# Plot trajectories of subjects in two dimensions

`plotPath()` joins each subject's samples in time order in a
two-dimensional space: two components of an ordination, or two variables
such as features or colData columns.

## Usage

``` r
plotPath(x, ...)

# S4 method for class 'SummarizedExperiment'
plotPath(
  x,
  time.col,
  group = NULL,
  dimred = NULL,
  ncomponents = 2,
  x.var = NULL,
  y.var = NULL,
  assay.type = "counts",
  origin = c("none", "baseline"),
  arrows = TRUE,
  colour.by = NULL,
  group.by = NULL,
  centroids = FALSE,
  field = NULL,
  lag = NULL,
  animate = FALSE,
  max.gap = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotPath(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_path`](https://ggplot2.tidyverse.org/reference/geom_path.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits. Time differences, as in `max.gap`, are in days for `Date` and
  hours for `POSIXct` times.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects. If `NULL`, all samples form one path. (Default:
  `NULL`)

- dimred:

  `Character scalar`. Specifies the `reducedDim(x)` to plot. (Default:
  `NULL`)

- ncomponents:

  `Numeric scalar` or `vector`. The two components of `dimred` to plot:
  `2` for the first two, or two components by index or name. (Default:
  `2`)

- x.var, y.var:

  `Character scalar`. Variables to plot instead of `dimred`: a column of
  `colData(x)` or a feature of `assay.type`. (Default: `NULL`)

- assay.type:

  `Character scalar`. Specifies the assay that holds features named in
  `x.var` or `y.var`. (Default: `"counts"`)

- origin:

  `Character scalar`. Whether paths start where the first sample lies
  (`"none"`) or at zero, by subtracting each subject's first sample
  (`"baseline"`). (Default: `"none"`)

- arrows:

  `Logical scalar`. Whether paths end in an arrow pointing forward in
  time. (Default: `TRUE`)

- colour.by:

  `Character scalar`. Specifies the column from `colData(x)` that
  colours the paths. (Default: `NULL`)

- group.by:

  `Character scalar`. Specifies the column from `colData(x)` that groups
  subjects for `centroids`. (Default: `NULL`)

- centroids:

  `Logical scalar`. Draws the mean path of each `group.by` group, or of
  all subjects: the mean position at each time. Best with times shared
  by the subjects, such as visits. (Default: `FALSE`)

- field:

  `Integer scalar`. Draws the mean movement to the next sample from each
  cell of a `field` by `field` grid as an arrow. If `NULL`, draws no
  field. (Default: `NULL`)

- lag:

  `Integer scalar`. Instead of paths, draws the value of `x.var`, or the
  first component of `dimred`, against its value `lag` samples later in
  the same path, with the identity line. (Default: `NULL`)

- animate:

  `Logical scalar`. Reveals the paths over time with
  [`gganimate::transition_reveal()`](https://gganimate.com/reference/transition_reveal.html);
  requires gganimate. (Default: `FALSE`)

- max.gap:

  `Numeric scalar`. Breaks a path where consecutive samples of a subject
  are further apart than this, in `time.col` units. For visits, `NULL`
  breaks paths at visits without a sample; otherwise `NULL` never breaks
  paths. The movement field and lag plots pair samples within paths
  only. (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Give either `dimred` or both `x.var` and `y.var`; with `lag`, `x.var`
alone is enough. Points near the identity line of a lag plot change
little over `lag` samples. Ordinations and lag plots have equal axis
scales.

## Examples

``` r
data("hitchip1006", package = "miaTime")
tse <- hitchip1006
# Some profiles are stored twice; keep the first copy
tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
# Subjects sampled five times
tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]

# Paths of subjects in a PCoA of Bray-Curtis dissimilarities
tse <- mia::transformAssay(tse, method = "relabundance")
tse <- mia::addMDS(tse, assay.type = "relabundance", method = "bray",
    name = "PCoA")
plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
    colour.by = "subject")


# The same paths, all starting at zero, to compare their directions
plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
    origin = "baseline", colour.by = "subject")


# Two genera against each other
plotPath(tse, time.col = "time", group = "subject",
    x.var = "Bacteroides vulgatus et rel.",
    y.var = "Prevotella melaninogenica et rel.",
    assay.type = "relabundance", colour.by = "subject")


# Mean movement over a 6 x 6 grid of the PCoA
plotPath(tse, time.col = "time", group = "subject", dimred = "PCoA",
    field = 6, arrows = FALSE)


# Mean path per visit and sex
tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
plotPath(tse, time.col = "visit", group = "subject", dimred = "PCoA",
    group.by = "sex", centroids = TRUE)


# Diversity against diversity at the next sample
plotPath(tse, time.col = "time", group = "subject", x.var = "diversity",
    lag = 1)
```
