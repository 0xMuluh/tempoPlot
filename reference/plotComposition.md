# Plot composition over time

`plotComposition()` stacks the relative abundances of the most abundant
features over time, for each subject or averaged over groups of
subjects.

## Usage

``` r
plotComposition(x, ...)

# S4 method for class 'SummarizedExperiment'
plotComposition(
  x,
  time.col,
  group = NULL,
  assay.type = "counts",
  features = NULL,
  top = 5,
  as.relative = TRUE,
  layout = c("bar", "area", "stream"),
  level = c("subject", "group"),
  group.by = NULL,
  max.gap = NULL,
  ...
)

# S4 method for class 'MultiAssayExperiment'
plotComposition(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_col`](https://ggplot2.tidyverse.org/reference/geom_bar.html) or
  [`geom_ribbon`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling times: numeric, `Date`, `POSIXct` or a factor of
  visits.

- group:

  `Character scalar`. Specifies the column from `colData(x)` that
  identifies subjects. Each subject gets its own panel. (Default:
  `NULL`)

- assay.type:

  `Character scalar`. Specifies the assay of abundances. To plot
  taxonomic ranks, agglomerate `x` first with
  [`mia::agglomerateByRank()`](https://microbiome.github.io/mia/reference/agglomerate-methods.html).
  (Default: `"counts"`)

- features:

  `Character vector`. Features to plot. If `NULL`, plots the `top`
  features with the largest mean plotted value: proportion, or assay
  value with `as.relative = FALSE`. The rest, and a feature already
  named "Other", are summed into "Other". (Default: `NULL`)

- top:

  `Integer scalar`. Number of features to plot when `features` is
  `NULL`. (Default: `5`)

- as.relative:

  `Logical scalar`. Turns the assay into proportions of each sample
  (`TRUE`) or stacks its values as they are (`FALSE`), for absolute
  abundances such as relative abundances times a qPCR or flow cytometry
  load. (Default: `TRUE`)

- layout:

  `Character scalar`. Draws stacked bars (`"bar"`), stacked areas
  (`"area"`) or areas stacked around zero (`"stream"`), which shows
  changes of the total with `as.relative = FALSE`. (Default: `"bar"`)

- level:

  `Character scalar`. Plots each subject (`"subject"`) or the mean
  composition of each `group.by` group at each time (`"group"`).
  (Default: `"subject"`)

- group.by:

  `Character scalar`. Specifies the column from `colData(x)` that groups
  subjects when `level = "group"`. If `NULL`, all subjects form one
  group. (Default: `NULL`)

- max.gap:

  `Numeric scalar`. Breaks areas and streams where consecutive times of
  a panel are further apart than this, in `time.col` units. For visits,
  `NULL` breaks them at visits without a sample; otherwise `NULL` never
  breaks them. A time between two breaks draws no area; bars show it.
  (Default: `NULL`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Examples

``` r
# A defined community of 16 species in three replicate bioreactors
data("minimalgut", package = "miaTime")
tse <- minimalgut
plotComposition(tse, time.col = "Time_hr_num",
    group = "StudyIdentifier", top = 8, layout = "area")


# Biomass of each species, relative abundance times optical density,
# as a stream
tse <- tse[, !is.na(tse$OpticalDensity)]
tse <- mia::transformAssay(tse, method = "relabundance")
SummarizedExperiment::assay(tse, "biomass") <-
    t(t(SummarizedExperiment::assay(tse, "relabundance")) *
    tse$OpticalDensity)
plotComposition(tse, time.col = "Time_hr_num",
    group = "StudyIdentifier", assay.type = "biomass", top = 8,
    as.relative = FALSE, layout = "stream")


# Phyla of five subjects over their visits
data("hitchip1006", package = "miaTime")
tse <- hitchip1006
# Some profiles are stored twice; keep the first copy
tse <- tse[, !duplicated(t(SummarizedExperiment::assay(tse, "counts")))]
tse <- tse[, tse$subject %in% names(which(table(tse$subject) == 5))]
tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
tse <- mia::agglomerateByRank(tse, rank = "Phylum")
plotComposition(tse, time.col = "visit", group = "subject")


# Mean phylum composition by sex at each visit
plotComposition(tse, time.col = "visit", group = "subject",
    level = "group", group.by = "sex")
```
