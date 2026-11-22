# Shared look of the plot builders: one theme and four colour scales.
# Not exported. Users override with + theme() or + scale_*(); ggplot2 then
# reports that the scale is replaced.

# Missing values in every scale: a mid grey, distinct from the white panel
# and from the light grey of zero on the diverging scale.
.na_colour <- "grey50"

# Okabe-Ito without yellow and black, in a fixed order. Every prefix passes
# the colour-vision check on a white panel (all pairs; worst deutan
# Delta E 7.6, first three 11.4).
.okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
    "#D55E00")

# White panel, light grey major grid, legend below the panel so that long
# time axes keep the full width.
#' @importFrom ggplot2 theme_bw theme element_blank
.theme_tempoPlot <- function() {
    theme_bw() +
        theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

# n categorical colours: Okabe-Ito in order. With more than six levels the
# whole scale switches to viridis, so that no colour repeats.
#' @importFrom grDevices hcl.colors
.palette_categorical <- function(n) {
    if (n <= length(.okabe_ito)) .okabe_ito[seq_len(n)] else
        hcl.colors(n, "viridis")
}

#' @importFrom ggplot2 discrete_scale
.scale_categorical <- function(aesthetics, ...) {
    discrete_scale(aesthetics, palette = .palette_categorical,
        na.value = .na_colour, ...)
}

#' @importFrom ggplot2 scale_colour_viridis_c
.scale_sequential <- function(aesthetics, ...) {
    scale_colour_viridis_c(..., na.value = .na_colour,
        aesthetics = aesthetics)
}

# Centred on 0: blue below, red above (the direction of miaViz), light grey
# at 0.
#' @importFrom ggplot2 scale_colour_gradient2
.scale_diverging <- function(aesthetics, ...) {
    scale_colour_gradient2(..., low = "#2166AC", mid = "grey90",
        high = "#B2182B", midpoint = 0, na.value = .na_colour,
        aesthetics = aesthetics)
}

# Scale for the column col of df (colour.by, may be NULL): sequential for
# numbers, else categorical.
.scale_by <- function(df, col, aesthetics, ...) {
    if (!is.null(col) && is.numeric(df[[col]])) {
        .scale_sequential(aesthetics, ...)
    } else {
        .scale_categorical(aesthetics, ...)
    }
}
