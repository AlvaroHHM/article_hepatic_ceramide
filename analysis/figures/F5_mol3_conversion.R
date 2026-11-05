# ============================================================
# F5_mol3_conversion.R
# Figure 5: MOL3_EAE sphingolipid genes — fraction vs intensity.
# Panel A: % positive cells (mean ± SD per mouse)
# Panel B: mean CPM in positive cells (mean ± SD per mouse)
# Error bars: SD per mouse (from 36c)
# Asterisks: per-mouse glmer binomial (from 36d)
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

if (dir.exists("article_hepatic_ceramide")) {
  REPO <- "article_hepatic_ceramide"
} else if (basename(getwd()) == "article_hepatic_ceramide") {
  REPO <- "."
} else stop("Run from ~/HARDVARD/recount3_analysis/ or from repo root")

PARENT <- if (REPO == ".") ".." else "."
source(file.path(REPO, "analysis/figures/00_figure_style.R"))
message("[style] loaded")

# ---- Load 36c plot data (means + SD) ----
cand_36c <- c(
  file.path(REPO, "results/snrnaseq/36c_positivity_plot_data.csv"),
  file.path(PARENT, "results/snrnaseq/36c_positivity_plot_data.csv")
)
path_36c <- cand_36c[file.exists(cand_36c)][1]
if (is.na(path_36c)) stop("36c_positivity_plot_data.csv not found")
d <- read.csv(path_36c, check.names = FALSE)
message("[data 36c] ", path_36c)

# ---- Load 36d glmer for asterisks ----
cand_36d <- c(
  file.path(REPO, "results/snrnaseq/36d_glmer_per_mouse.csv"),
  file.path(PARENT, "results/snrnaseq/36d_glmer_per_mouse.csv")
)
path_36d <- cand_36d[file.exists(cand_36d)][1]
if (is.na(path_36d)) stop("36d_glmer_per_mouse.csv not found")
d_glmer <- read.csv(path_36d, check.names = FALSE)
message("[data 36d] ", path_36d)

# ---- Merge + asterisks from glmer ----
genes <- c("B4galt6", "Ugcg", "Smpd3", "Sgpl1", "S1pr1")
d <- d %>%
  filter(gene %in% genes) %>%
  left_join(d_glmer %>% dplyr::select(gene, glmer_p_adj), by = "gene") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Early", "Peak", "Late")),
    ast = dplyr::case_when(
      is.na(glmer_p_adj)  ~ "ns",
      glmer_p_adj < 0.001 ~ "***",
      glmer_p_adj < 0.01  ~ "**",
      glmer_p_adj < 0.05  ~ "*",
      glmer_p_adj < 0.10  ~ "#",
      TRUE                ~ "ns"
    )
  )
message("[ast] per gene: ",
        paste(unique(paste0(d$gene, "=", d$ast)), collapse = ", "))

annot <- d %>% filter(timepoint == "Peak")

# ---- Panel A: fraction (%) ----
pA <- ggplot(d, aes(x = timepoint, y = mean_pct, group = gene, color = gene)) +
  geom_errorbar(aes(ymin = pmax(mean_pct - sd_pct, 0),
                    ymax = mean_pct + sd_pct),
                width = 0.15, linewidth = 0.4, alpha = 0.7) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  geom_text(data = annot,
            aes(x = timepoint, y = mean_pct + sd_pct, label = ast),
            vjust = -0.8, size = 3.5, fontface = "bold",
            show.legend = FALSE) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(
    "B4galt6" = "#EE0000", "Ugcg" = "#3B4992", "Smpd3" = "#008B45",
    "Sgpl1"   = "#631879", "S1pr1" = "#F0A202"
  )) +
  theme_pub(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "% positive cells")

# ---- Panel B: intensity (CPM) ----
pB <- ggplot(d, aes(x = timepoint, y = mean_cpm, group = gene, color = gene)) +
  geom_errorbar(aes(ymin = pmax(mean_cpm - sd_cpm, 0),
                    ymax = mean_cpm + sd_cpm),
                width = 0.15, linewidth = 0.4, alpha = 0.7) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  geom_text(data = annot,
            aes(x = timepoint, y = mean_cpm + sd_cpm, label = ast),
            vjust = -0.8, size = 3.5, fontface = "bold",
            show.legend = FALSE) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(
    "B4galt6" = "#EE0000", "Ugcg" = "#3B4992", "Smpd3" = "#008B45",
    "Sgpl1"   = "#631879", "S1pr1" = "#F0A202"
  )) +
  theme_pub(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "mean CPM (positive cells)")

library(cowplot)
final <- plot_grid(pA, pB, ncol = 1, labels = c("A", "B"),
                   label_size = 11, align = "v")

out_dir <- file.path(REPO, "figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "F5_mol3_conversion.pdf")
ggsave(out_pdf, final, width = 10, height = 5.5, dpi = 300)
message("[OK] ", out_pdf)
