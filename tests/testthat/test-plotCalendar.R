library(SummarizedExperiment)
library(ggplot2)

# 2024-01-01 is a Monday.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 5L)),
    colData = DataFrame(d = as.Date(c("2024-01-01", "2024-01-01",
        "2024-01-03", "2024-01-08", "2024-12-30")), n = 1:5))

test_that("plotCalendar checks its arguments", {
    expect_error(plotCalendar(toy, "n"), "Date")
    expect_error(plotCalendar(toy, "none"), "'time.col'")
    expect_error(plotCalendar(toy, "d", unit = "month"))
})

test_that("plotCalendar counts samples per day", {
    df <- plotCalendar(toy, "d")$data
    expect_identical(nrow(df), 4L)
    jan1 <- df[df$.x == as.Date("2024-01-01") & df$.y == "Mon", ]
    expect_identical(jan1$n, 2L)
    expect_identical(as.character(df$.y[df$n == 1L & df$.x ==
        as.Date("2024-01-01")]), "Wed")
})

test_that("plotCalendar counts samples per ISO week", {
    df <- plotCalendar(toy, "d", unit = "week")$data
    # Week 1 of 2024 has three samples, week 2 one; 2024-12-30 is in
    # ISO week 1 of 2025.
    expect_identical(df$n[df$.y == "2024" & df$.x == 1L], 3L)
    expect_identical(df$n[df$.y == "2024" & df$.x == 2L], 1L)
    expect_identical(df$n[df$.y == "2025" & df$.x == 1L], 1L)
})

test_that("plotCalendar fills tiles on the sequential scale", {
    p <- plotCalendar(toy, "d")
    ld <- layer_data(p)
    # Fewest samples darkest, most brightest (viridis).
    expect_identical(unique(ld$fill[p$data$n == 1L]), "#440154")
    expect_identical(unique(ld$fill[p$data$n == 2L]), "#FDE725")
    expect_identical(p$theme$legend.position, "bottom")
})

test_that("plotCalendar passes ... to geom_tile", {
    p <- plotCalendar(toy, "d", colour = "black")
    expect_true(all(layer_data(p)$colour == "black"))
})
