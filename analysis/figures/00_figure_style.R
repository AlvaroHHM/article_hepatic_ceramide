# ============================================================
# 00_figure_style.R
# Shared publication style for all manuscript figures.
#   - theme_cowplot (from cowplot) as base
#   - AAAS palette (ggsci::pal_aaas) hardcoded for portability
#   - No titles inside plots (captions in manuscript)
#   - Panel labels via cowplot::plot_grid(labels = c("A)", ...))
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

# AAAS palette (from ggsci::pal_aaas)
AAAS <- list(
  blue    = "#3B4992",
  red     = "#EE0000",
  green   = "#008B45",
  purple  = "#631879",
  teal    = "#008280",
  darkred = "#BB0021",
  grey    = "#808180",
  black   = "#1B1919"
)

# Publication-grade theme: clean, no grid, axis lines only
theme_pub <- function(base_size = 9) {
  theme_cowplot(font_size = base_size, line_size = 0.4) +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      axis.line         = element_line(color = "black", linewidth = 0.3),
      axis.ticks        = element_line(color = "black", linewidth = 0.3),
      axis.text         = element_text(color = "black", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size),
      legend.position   = "right",
      legend.key.size   = unit(0.35, "cm"),
      legend.text       = element_text(size = base_size - 1),
      legend.title      = element_text(size = base_size, face = "bold"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA)
    )
}

# AAAS blue-white-red gradient (low = blue, high = red)
scale_fill_aaas_gradient <- function(midpoint = 0, name = "log2FC",
                                       limits = NULL) {
  scale_fill_gradient2(
    low      = AAAS$blue,
    mid      = "#FFFFFF",
    high     = AAAS$red,
    midpoint = midpoint,
    name     = name,
    limits   = limits,
    na.value = "grey90"
  )
}

# Same gradient but for colors (lines/points)
scale_color_aaas_gradient <- function(midpoint = 0, name = "log2FC") {
  scale_color_gradient2(
    low      = AAAS$blue,
    mid      = "#FFFFFF",
    high     = AAAS$red,
    midpoint = midpoint,
    name     = name
  )
}

message("[style] theme_pub, AAAS palette, aaas_gradient loaded")
