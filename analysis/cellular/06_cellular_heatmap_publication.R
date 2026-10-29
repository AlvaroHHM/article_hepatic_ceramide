if (dir.exists("article_hepatic_ceramide")) {
  source("article_hepatic_ceramide/analysis/figures/00_figure_style.R")
} else if (file.exists("analysis/figures/00_figure_style.R")) {
  source("analysis/figures/00_figure_style.R")
}
# ============================================================
# 06_cellular_heatmap_publication.R
# F3: Sphingolipid genes across CNS cell types (basal)
# - No title inside plot (goes in manuscript caption)
# - AAAS palette
# - Category labels without grey background boxes
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble)
})

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

long_path <- "results/cellular/27_cellular_comparison_long.csv"
if (!file.exists(long_path)) {
  long_path <- "results/27_cellular_comparison_long.csv"
}
if (!file.exists(long_path)) {
  stop("Falta 27_cellular_comparison_long.csv")
}
long_df <- read.csv(long_path, stringsAsFactors = FALSE)

cells_expected <- c("Astrocytes", "Microglia", "Glia_NEG", "Neurons_POS")

genes_ok <- long_df %>%
  filter(!is.na(mean_Healthy), cell %in% cells_expected) %>%
  count(gene, name = "n_cells") %>%
  filter(n_cells == length(cells_expected)) %>%
  pull(gene)

if (length(genes_ok) == 0) stop("No gene has data in all 4 cells")

gene_order <- c("Sgms1","Sgms2","Smpd3","Ugcg","B4galt5","B4galt6",
                "Asah1","Sgpl1","S1pr1","S1pr3","Pla2g4a")
gene_order <- gene_order[gene_order %in% genes_ok]
extras <- setdiff(genes_ok, gene_order)
if (length(extras) > 0) gene_order <- c(gene_order, sort(extras))

cell_order <- cells_expected[cells_expected %in% long_df$cell]

mat <- long_df %>%
  filter(gene %in% gene_order, cell %in% cell_order, !is.na(mean_Healthy)) %>%
  select(gene, cell, mean_Healthy) %>%
  mutate(gene = factor(gene, levels = rev(gene_order)),
         cell = factor(cell, levels = cell_order))

cat_map <- c(
  "Sgms1" = "SMS (Cer->SM)",  "Sgms2" = "SMS (Cer->SM)",
  "Smpd3" = "SMPD (SM->Cer)",
  "Ugcg" = "LacCer synthesis", "B4galt5" = "LacCer synthesis",
  "B4galt6" = "LacCer synthesis",
  "Asah1" = "Ceramide degradation", "Sgpl1" = "Ceramide degradation",
  "S1pr1" = "S1P receptor", "S1pr3" = "S1P receptor",
  "Pla2g4a" = "cPLA2 effector"
)
mat$category <- cat_map[as.character(mat$gene)]
mat$category[is.na(mat$category)] <- "Other"

# ---- HEATMAP 1: absolute log2(RPKM+1), AAAS ----
mat$log_val <- log2(mat$mean_Healthy + 1)

p_abs <- ggplot(mat, aes(x = cell, y = gene, fill = log_val)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient2(low      = "#3B4992",
                       mid      = "#FFFFFF",
                       high     = "#EE0000",
                       midpoint = median(mat$log_val, na.rm = TRUE),
                       name     = "log2(RPKM+1)") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  theme_pub(base_size = 10) +
  theme(
    axis.text.x    = element_text(angle = 30, hjust = 1, size = 9, face = "bold"),
    axis.text.y    = element_text(size = 9),
    strip.background = element_blank(),
    strip.text.y   = element_text(angle = 0, size = 7, face = "bold", hjust = 0),
    panel.grid     = element_blank()
  ) +
  labs(x = NULL, y = NULL)

ggsave("results/cellular/27b_cellular_heatmap_absolute.pdf",
       p_abs, width = 7.5, height = 6.5)
cat("[OK] 27b_cellular_heatmap_absolute.pdf\n")

# ---- HEATMAP 2: z-score within column, AAAS ----
mat_z <- mat %>%
  group_by(cell) %>%
  mutate(z = if (n() > 1 && sd(mean_Healthy, na.rm = TRUE) > 0)
            (mean_Healthy - mean(mean_Healthy, na.rm = TRUE)) /
              sd(mean_Healthy, na.rm = TRUE) else 0) %>%
  ungroup()
mat_z$z_clip <- pmax(pmin(mat_z$z, 3), -3)

p_z <- ggplot(mat_z, aes(x = cell, y = gene, fill = z_clip)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient2(low      = "#3B4992",
                       mid      = "#FFFFFF",
                       high     = "#EE0000",
                       midpoint = 0, limits = c(-3, 3),
                       name     = "z-score\n(within column)") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  theme_pub(base_size = 10) +
  theme(
    axis.text.x    = element_text(angle = 30, hjust = 1, size = 9, face = "bold"),
    axis.text.y    = element_text(size = 9),
    strip.background = element_blank(),
    strip.text.y   = element_text(angle = 0, size = 7, face = "bold", hjust = 0),
    panel.grid     = element_blank()
  ) +
  labs(x = NULL, y = NULL)

ggsave("results/cellular/27b_cellular_heatmap_zscore.pdf",
       p_z, width = 7.5, height = 6.5)
cat("[OK] 27b_cellular_heatmap_zscore.pdf\n")
