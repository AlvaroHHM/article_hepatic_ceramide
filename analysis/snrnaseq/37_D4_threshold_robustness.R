suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
setwd("~/HARDVARD/recount3_analysis")

# Try multiple locations
ckpt_candidates <- c(
  "results/35b_checkpoint",
  "article_hepatic_ceramide/results/35b_checkpoint",
  "../results/35b_checkpoint"
)
ckpt_dir <- ckpt_candidates[dir.exists(ckpt_candidates)][1]
if (is.na(ckpt_dir)) {
  cat("[FAIL] No encontré 35b_checkpoint en:\n")
  for (p in ckpt_candidates) cat("  ", p, "(exists:", dir.exists(p), ")\n")
  quit(status = 1)
}

files <- list.files(ckpt_dir, pattern = "\\.rds$", full.names = TRUE)
cat("[OK] Checkpoint dir:", ckpt_dir, "| files:", length(files), "\n")

cells <- bind_rows(lapply(files, readRDS))
cat("[OK] Total cells:", nrow(cells), "\n")
cat("[OK] Conditions:", paste(unique(cells$condition), collapse = ", "), "\n")
cat("[OK] Timepoints:", paste(unique(cells$timepoint), collapse = ", "), "\n")
cat("[OK] Genes available:", paste(intersect(c("B4galt6","Ugcg","Sgpl1","Smpd3","S1pr1"), colnames(cells)), collapse = ", "), "\n")

# Per-mouse fraction at each threshold, then mean across mice
genes <- c("B4galt6", "Ugcg", "Sgpl1", "Smpd3", "S1pr1")
thresholds <- c(0, 0.5, 1, 2)
results <- data.frame()

for (thr in thresholds) {
  for (g in genes) {
    if (!g %in% colnames(cells)) next

    per_mouse <- cells %>%
      filter(condition == "EAE", timepoint %in% c("Early", "Peak")) %>%
      group_by(orig.ident, timepoint) %>%
      summarise(
        frac = mean(.data[[g]] > thr, na.rm = TRUE),
        n_cells = n(),
        .groups = "drop"
      ) %>%
      group_by(timepoint) %>%
      summarise(
        mean_frac = mean(frac, na.rm = TRUE),
        n_mice = n(),
        .groups = "drop"
      )

    if (nrow(per_mouse) == 2 &&
        all(c("Early", "Peak") %in% per_mouse$timepoint)) {
      fE <- per_mouse$mean_frac[per_mouse$timepoint == "Early"]
      fP <- per_mouse$mean_frac[per_mouse$timepoint == "Peak"]
      results <- rbind(results, data.frame(
        gene = g,
        threshold = thr,
        frac_Early = round(fE, 4),
        frac_Peak = round(fP, 4),
        ratio_peak_early = round(fP / fE, 3),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n=== Threshold robustness of conversion ===\n")
print(results, row.names = FALSE)

# Interpretación automática
cat("\n=== INTERPRETATION ===\n")
for (g in genes) {
  sub <- results[results$gene == g, ]
  if (nrow(sub) == 0) next
  robust <- all(sub$ratio_peak_early > 1.5, na.rm = TRUE)
  cat(sprintf("  %-8s ratio range [%.2f, %.2f] across thresholds → %s\n",
              g, min(sub$ratio_peak_early), max(sub$ratio_peak_early),
              if (robust) "ROBUST (>1.5 at all thr)" else "SENSITIVE to threshold"))
}

dir.create("results/snrnaseq", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "results/snrnaseq/37_D4_threshold_robustness.csv", row.names = FALSE)
cat("\n[OK] results/snrnaseq/37_D4_threshold_robustness.csv\n")
