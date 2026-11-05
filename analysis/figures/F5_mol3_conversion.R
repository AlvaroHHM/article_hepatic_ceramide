# ============================================================
# F5_mol3_conversion.R
# Figure 5: MOL3_EAE sphingolipid genes — fraction vs intensity.
#
# Panel A: % positive cells. Jitter per mouse (n=5-6) + mean line +
#          SEM error bars + asterisks (per-mouse glmer, only p_adj<0.05).
# Panel B: mean CPM in positive cells. Jitter per mouse + mean line +
#          SEM error bars. NO asterisks (intensity was not tested; the
#          glmer binomial tests fraction positivity only).
#
# Inputs:
#   results/snrnaseq/36c_positivity_plot_data.csv    (means + SD)
#   results/snrnaseq/36c_positivity_per_mouse.csv    (per-mouse, jitter)
#   results/snrnaseq/36d_glmer_per_mouse.csv         (asterisks)
# Output: figures/main/F5_mol3_conversion.pdf
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
d_agg <- read.csv(path_36c, check.names = FALSE)
message("[data 36c agg] ", path_36c)

# ---- Load 36c per-mouse data (for jitter) ----
cand_mouse <- c(
  file.path(REPO, "results/snrnaseq/36c_positivity_per_mouse.csv"),
  file.path(PARENT, "results/snrnaseq/36c_positivity_per_mouse.csv")
)
path_mouse <- cand_mouse[file.exists(cand_mouse)][1]
if (is.na(path_mouse)) stop("36c_positivity_per_mouse.csv not found")
d_mouse <- read.csv(path_mouse, check.names = FALSE)
message("[data 36c mouse] ", path_mouse)

# ---- Load 36d glmer for asterisks ----
cand_36d <- c(
  file.path(REPO, "results/snrnaseq/36d_glmer_per_mouse.csv"),
  file.path(PARENT, "results/snrnaseq/36d_glmer_per_mouse.csv")
)
path_36d <- cand_36d[file.exists(cand_36d)][1]
if (is.na(path_36d)) stop("36d_glmer_per_mouse.csv not found")
d_glmer <- read.csv(path_36d, check.names = FALSE)
message("[data 36d] ", path_36d)

# ---- Merge + compute SEM and ast_panelA (only sig) ----
genes <- c("B4galt6", "Ugcg", "Smpd3", "Sgpl1", "S1pr1")

d_agg <- d_agg %>%
  filter(gene %in% genes) %>%
  left_join(d_glmer %>% dplyr::select(gene, glmer_p_adj), by = "gene") %>%
  mutate(
    timepoint  = factor(timepoint, levels = c("Early", "Peak", "Late")),
    sem_pct    = sd_pct / sqrt(n_mice),
    sem_cpm    = sd_cpm / sqrt(n_mice),
    ast_panelA = dplyr::case_when(
      is.na(glmer_p_adj)  ~ "",
      glmer_p_adj < 0.001 ~ "***",
      glmer_p_adj < 0.01  ~ "**",
      glmer_p_adj < 0.05  ~ "*",
      TRUE                ~ ""
    )
  )

d_mouse <- d_mouse %>%
  filter(gene %in% genes) %>%
  mutate(timepoint = factor(timepoint, levels = c("Early", "Peak", "Late")))

message("[ast sig-only] ",
        paste(unique(paste0(d_agg$gene, "=",
                            ifelse(d_agg$ast_panelA == "", "ns", d_agg$ast_panelA))),
              collapse = ", "))

annot <- d_agg %>% filter(timepoint == "Peak", ast_panelA != "")

# ---- Panel A: fraction (%) ----
pA <- ggplot(d_agg, aes(x = timepoint, y = mean_pct, group = gene, color = gene)) +
  geom_jitter(data = d_mouse,
              aes(x = timepoint, y = pct_positive),
              width = 0.10, size = 1.3, alpha = 0.45, color = "gray40",
              inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = pmax(mean_pct - sem_pct, 0),
                    ymax = mean_pct + sem_pct),
                width = 0.15, linewidth = 0.5, alpha = 0.85) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  geom_text(data = annot,
            aes(x = timepoint, y = mean_pct + sem_pct, label = ast_panelA),
            vjust = -0.8, size = 3.5, fontface = "bold",
            inherit.aes = FALSE, color = "black") +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(
    "B4galt6" = "#EE0000", "Ugcg" = "#3B4992", "Smpd3" = "#008B45",
    "Sgpl1"   = "#631879", "S1pr1" = "#F0A202"
  )) +
  theme_pub(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "% positive cells")

# ---- Panel B: intensity (CPM), NO asterisks ----
pB <- ggplot(d_agg, aes(x = timepoint, y = mean_cpm, group = gene, color = gene)) +
  geom_jitter(data = d_mouse,
              aes(x = timepoint, y = mean_pos_CPM),
              width = 0.10, size = 1.3, alpha = 0.45, color = "gray40",
              inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = pmax(mean_cpm - sem_cpm, 0),
                    ymax = mean_cpm + sem_cpm),
                width = 0.15, linewidth = 0.5, alpha = 0.85) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
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
ggsave(out_pdf, final, width = 10, height = 6.5, dpi = 300)
message("[OK] ", out_pdf)
