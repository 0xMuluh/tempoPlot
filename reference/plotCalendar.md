# Plot sample counts on a calendar

`plotCalendar()` counts the samples taken on each day or in each week
and draws them as calendar tiles: weeks by weekday, or years by week of
the year.

## Usage

``` r
plotCalendar(x, ...)

# S4 method for class 'SummarizedExperiment'
plotCalendar(x, time.col, unit = c("day", "week"), ...)

# S4 method for class 'MultiAssayExperiment'
plotCalendar(x, ..., experiment = 1L)
```

## Arguments

- x:

  a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- ...:

  additional arguments passed to
  [`geom_tile`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

- time.col:

  `Character scalar`. Specifies the column from `colData(x)` that
  contains sampling dates, of class `Date`, or times of class `POSIXct`,
  counted on their date in their time zone.

- unit:

  `Character scalar`. Counts samples per day, with weeks on the x axis
  and weekdays on the y axis (`"day"`), or per week, with the ISO week
  on the x axis and the year on the y axis (`"week"`). (Default:
  `"day"`)

- experiment:

  `Character scalar` or `integer scalar`. For a `MultiAssayExperiment`,
  the experiment to plot; its samples get the `colData` of `x`.
  (Default: `1L`)

## Value

A `ggplot` object.

## Details

Weeks start on Monday. Days and weeks without samples are left blank.

## Examples

``` r
# Daily sampling of 20 people; study days placed on an arbitrary
# calendar starting on 1 March, with no samples at weekends
data("temporalMicrobiome20", package = "miaTime")
tse <- temporalMicrobiome20
tse$date <- as.Date("2021-03-01") + tse$day - 1
tse <- tse[, !format(tse$date, "%u") %in% c("6", "7")]
plotCalendar(tse, time.col = "date")


# Samples per week of the year
plotCalendar(tse, time.col = "date", unit = "week")
```
