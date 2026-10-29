# ============================================================
# S2_MOL3_expression.R
# Supplementary S2: MOL3_EAE expression by timepoint
# Reconstructed (original 35b source was lost).
# Reads results/35b_MOL3_CPM_by_timepoint.csv and regenerates the plot
# with AAAS palette, no title (goes in caption).
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})

# Detect repo root
if (dir.exists("article_hepatic_ceramide")) {
  REPO <- "article_hepatic_ceramide"; PARENT <- "."
} else if (basename(getwd()) == "article_hepatic_ceramide") {
  REPO <- "."; PARENT <- ".."
} else stop("Run from ~/HARDVARD/recount3_analysis/ or repo root")

# Load style
style_path <- file.path(REPO, "analysis/figures/00_figure_style.R")
if (file.exists(style_path)) {
  source(style_path)
} else {
  AAAS <- list(blue="#3B4992", red="#EE0000", green="#008B45",
               purple="#631879", teal="#008280", grey="#808180")
  theme_pub <- function(base_size=9) {
    theme_cowplot(font_size=base_size) +
      theme(panel.grid=element_blank(),
            strip.background=element_blank())
  }
}
message("[style] loaded")

# Find data CSV
candidates <- c(
  file.path(REPO, "results/snrnaseq/35b_MOL3_CPM_by_timepoint.csv"),
  file.path(PARENT, "results/35b_MOL3_CPM_by_timepoint.csv"),
  file.path(REPO, "results/snrnaseq/35b_MOL3_early_vs_peak.csv"),
  file.path(PARENT, "results/35b_MOL3_early_vs_peak.csv")
)
path <- NA
for (c in candidates) if (file.exists(c)) { path <- c; break }
if (is.na(path)) stop("No 35b_MOL3 CSV found in expected locations")

message("[data] ", path)
d <- read.csv(path, check.names=FALSE)
message("[data] columns: ", paste(colnames(d), collapse=", "))
message("[data] rows: ", nrow(d))

# Detect format
if ("timepoint" %in% colnames(d) && "gene" %in% colnames(d)) {
  cpm_col <- intersect(c("mean_CPM","CPM_mean","CPM","mean","value"), colnames(d))[1]
  if (is.na(cpm_col)) stop("No CPM column. Cols: ", paste(colnames(d), collapse=", "))
  message("[format] long (using ", cpm_col, ")")
  plot_df <- d %>%
    select(gene, timepoint, value = all_of(cpm_col)) %>%
    filter(!is.na(value)) %>%
    mutate(timepoint = factor(timepoint, levels=c("Early","Peak","Late")))
} else {
  tp_cols <- intersect(c("Early","Peak","Late","early","peak","late"), colnames(d))
  if (length(tp_cols) < 2) stop("Cannot detect format. Cols: ", paste(colnames(d), collapse=", "))
  message("[format] wide (timepoints: ", paste(tp_cols, collapse=", "), ")")
  plot_df <- d %>%
    select(gene, all_of(tp_cols)) %>%
    pivot_longer(-gene, names_to="timepoint", values_to="value") %>%
    mutate(timepoint = factor(timepoint, levels=c("Early","Peak","Late")))
}

# Pick 8 genes with highest variance
genes_to_show <- plot_df %>%
  group_by(gene) %>%
  summarise(v = sd(value, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(v)) %>%
  head(8) %>%
  pull(gene)

message("[genes] ", paste(genes_to_show, collapse=", "))
plot_df <- plot_df %>% filter(gene %in% genes_to_show)

# AAAS color map
gene_levels <- sort(unique(plot_df$gene))
palette <- c(AAAS$red, AAAS$blue, AAAS$green, AAAS$purple, AAAS$teal,
             "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#999999")
color_map <- setNames(palette[seq_along(gene_levels)], gene_levels)

p <- ggplot(plot_df, aes(x=timepoint, y=value, color=gene, group=gene)) +
  geom_line(linewidth=0.9) +
  geom_point(size=2.5) +
  scale_color_manual(values=color_map, name=NULL) +
  facet_wrap(~gene, scales="free_y", nrow=2) +
  theme_pub(base_size=10) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face="bold", size=10),
    panel.grid       = element_blank(),
    legend.position  = "none"
  ) +
  labs(x=NULL, y="CPM (mean per timepoint)")

out_pdf <- file.path(REPO, "figures/supplementary/S2_MOL3_expression.pdf")
ggsave(out_pdf, p, width=10, height=5, dpi=300)
message("[OK] ", out_pdf)
