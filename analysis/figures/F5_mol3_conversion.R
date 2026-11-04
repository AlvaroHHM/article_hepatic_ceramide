# ============================================================
# F5_mol3_conversion.R
# Figure 5: MOL3_EAE B4galt6+ cells — fraction vs intensity.
# - Fraction (%) = pct_positive  -> panel A
# - Intensity (CPM) = mean_pos_CPM -> panel B
#
# Data: results/snrnaseq/36_positivity.csv (or ../results/36_positivity.csv)
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

# Locate input
candidates <- c(
  file.path(REPO, "results/snrnaseq/36_positivity.csv"),
  file.path(PARENT, "results/36_positivity.csv")
)
path <- NA
for (p in candidates) if (file.exists(p)) { path <- p; break }
if (is.na(path)) stop("36_positivity.csv not found")

message("[data] ", path)
d <- read.csv(path, check.names = FALSE)
message("[data] columns: ", paste(colnames(d), collapse=", "))

need <- c("gene", "timepoint", "pct_positive", "mean_pos_CPM")
miss <- setdiff(need, colnames(d))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse=", "))

genes <- c("B4galt6", "Ugcg", "Smpd3", "Sgpl1", "S1pr1")
d <- d %>%
  filter(gene %in% genes) %>%
  mutate(timepoint = factor(timepoint, levels = c("Early", "Peak", "Late")))

# Panel A: fraction (%)
pA <- ggplot(d, aes(x = timepoint, y = pct_positive,
                    group = gene, color = gene)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_color_manual(values = AAAS[c("red","blue","green","purple","orange")]) +
  theme_pub(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "% positive cells")

# Panel B: intensity (mean CPM in positive cells)
pB <- ggplot(d, aes(x = timepoint, y = mean_pos_CPM,
                    group = gene, color = gene)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_color_manual(values = AAAS[c("red","blue","green","purple","orange")]) +
  theme_pub(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "mean CPM (positive cells)")

# Combine vertically
library(cowplot)
final <- plot_grid(pA, pB, ncol = 1, labels = c("A", "B"),
                   label_size = 11, align = "v")

out_dir <- file.path(REPO, "figures/main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "F5_mol3_conversion.pdf")
ggsave(out_pdf, final, width = 10, height = 5.5, dpi = 300)
message("[OK] ", out_pdf)
