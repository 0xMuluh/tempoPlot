library(ggplot2)

df <- data.frame(x = 1:4, y = 1:4, g = c("a", "b", "a", NA),
    v = c(-1, 0, 1, NA))

test_that(".theme_tempoPlot has a white panel and a light major grid", {
    th <- .theme_tempoPlot()
    expect_identical(calc_element("panel.background", th)$fill, "white")
    expect_identical(col2rgb(calc_element("panel.grid.major", th)$colour),
        col2rgb("grey92"))
    expect_true(inherits(calc_element("panel.grid.minor", th),
        "element_blank"))
    expect_identical(th$legend.position, "bottom")
    th <- th + theme(legend.position = "right")
    expect_identical(th$legend.position, "right")
})

test_that(".palette_categorical is Okabe-Ito in order, then viridis", {
    expect_identical(.palette_categorical(2), c("#0072B2", "#E69F00"))
    expect_identical(.palette_categorical(6), .okabe_ito)
    seven <- .palette_categorical(7)
    expect_identical(seven, grDevices::hcl.colors(7, "viridis"))
    expect_length(unique(seven), 7L)
})

test_that(".scale_categorical colours levels in order and NA grey", {
    p <- ggplot(df, aes(x, y, colour = g)) + geom_point() +
        .scale_categorical("colour")
    expect_identical(layer_data(p)$colour,
        c("#0072B2", "#E69F00", "#0072B2", "grey50"))
    p <- ggplot(df, aes(x, y, fill = g)) + geom_col() +
        .scale_categorical("fill")
    expect_identical(layer_data(p)$fill[1:2], c("#0072B2", "#E69F00"))
    many <- data.frame(x = 1:8, g = letters[1:8])
    p <- ggplot(many, aes(x, x, colour = g)) + geom_point() +
        .scale_categorical("colour")
    expect_identical(layer_data(p)$colour, grDevices::hcl.colors(8,
        "viridis"))
})

test_that(".scale_categorical draws the grey level in neutral grey", {
    d <- data.frame(x = 1:3, g = factor(c("a", "Other", "b"),
        levels = c("a", "b", "Other")))
    p <- ggplot(d, aes(x, x, colour = g)) + geom_point() +
        .scale_categorical("colour", grey = "Other")
    expect_identical(layer_data(p)$colour, c("#0072B2", "grey50",
        "#E69F00"))
    key <- get_guide_data(p, "colour")
    expect_identical(key$.label, c("a", "b", "Other"))
    expect_identical(key$colour, c("#0072B2", "#E69F00", "grey50"))
    # Without the grey level, every level keeps its palette colour.
    p <- ggplot(d[-2L, ], aes(x, x, colour = g)) + geom_point() +
        .scale_categorical("colour", grey = "Other")
    expect_identical(layer_data(p)$colour, c("#0072B2", "#E69F00"))
})

test_that(".scale_highlight draws yes in the first colour, no grey", {
    d <- data.frame(x = 1:3, h = c(FALSE, TRUE, NA))
    p <- ggplot(d, aes(x, x, colour = h)) + geom_point() +
        .scale_highlight("colour")
    expect_identical(layer_data(p)$colour, c("grey80", "#0072B2",
        "grey50"))
    key <- get_guide_data(p, "colour")
    expect_identical(key$.label, c("yes", "no"))
    expect_identical(key$colour, c("#0072B2", "grey80"))
})

test_that(".scale_sequential is viridis with grey NA", {
    p <- ggplot(df, aes(x, y, fill = v)) + geom_tile() +
        .scale_sequential("fill")
    expect_true(inherits(p$scales$get_scales("fill"), "ScaleContinuous"))
    fill <- layer_data(p)$fill
    expect_identical(fill[1L], "#440154")
    expect_identical(fill[4L], "grey50")
})

test_that(".scale_diverging is blue below 0, light grey at 0, red above", {
    p <- ggplot(df, aes(x, y, fill = v)) + geom_tile() +
        .scale_diverging("fill")
    expect_identical(layer_data(p)$fill,
        c("#2166AC", "#E5E5E5", "#B2182B", "grey50"))
    # 0 stays light grey when the data are not symmetric around it.
    p <- ggplot(data.frame(v = c(-1, 0, 3)), aes(v, v, colour = v)) +
        geom_point() + .scale_diverging("colour")
    expect_identical(layer_data(p)$colour[2:3], c("#E5E5E5", "#B2182B"))
})

test_that(".scale_by picks sequential for numbers, else categorical", {
    expect_true(inherits(.scale_by(df, "v", "colour"), "ScaleContinuous"))
    expect_true(inherits(.scale_by(df, "g", "colour"), "ScaleDiscrete"))
    expect_true(inherits(.scale_by(df, NULL, "colour"), "ScaleDiscrete"))
    expect_identical(.scale_by(df, "g", "fill")$aesthetics, "fill")
})

test_that("users override the shared scales with + scale_*()", {
    p <- ggplot(df, aes(x, y, colour = g)) + geom_point() +
        .scale_categorical("colour")
    expect_message(p <- p + scale_colour_manual(values = c("red", "blue")),
        "already present")
    expect_identical(layer_data(p)$colour[1:2], c("red", "blue"))
    p <- ggplot(df, aes(x, y, fill = v)) + geom_tile() +
        .scale_diverging("fill")
    expect_message(p <- p + scale_fill_viridis_c(), "already present")
    expect_identical(layer_data(p)$fill[1L], "#440154")
})
