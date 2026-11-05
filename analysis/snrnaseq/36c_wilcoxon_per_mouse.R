# ============================================================
# 36c_wilcoxon_per_mouse.R
# Wilcoxon per-mouse (pseudobulk) for the 5 genes in Figure 5.
# Replaces the per-cell Wilcoxon in 35b_MOL3_wilcoxon.csv (which has
# pseudoreplication: n=1557 cells vs 3197 cells, not n=6 vs 5 mice).
#
# Inputs:
#   ../results/35b_checkpoint/*.rds (26 files, 13,741 MOL3_EAE cells)
# Outputs:
#   results/snrnaseq/36c_wilcoxon_per_mouse.csv
#   results/snrnaseq/36c_positivity_per_mouse.csv
#   results/snrnaseq/36c_positivity_plot_data.csv
#   logs/snrnaseq/36c_wilcoxon_per_mouse_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

dir.create("logs/snrnaseq", recursive = TRUE, showWarnings = FALSE)
sink("logs/snrnaseq/36c_wilcoxon_per_mouse_log.txt", split = TRUE)
on.exit(sink())

cat("=== 36c: Wilcoxon per-mouse for Figure 5 genes ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# ---- 1. Find checkpoints ----
ckpt_candidates <- c(
  "results/35b_checkpoint",
  "../results/35b_checkpoint"
)
ckpt_dir <- ckpt_candidates[dir.exists(ckpt_candidates)][1]
if (is.na(ckpt_dir)) stop("No checkpoint dir found. Looked in:\n",
                          paste("  -", ckpt_candidates, collapse = "\n"))
ckpts <- list.files(ckpt_dir, pattern = "\\.rds$", full.names = TRUE)
cat("[1/6] Checkpoints:", length(ckpts), "files in", ckpt_dir, "\n")
stopifnot(length(ckpts) > 0)

# ---- 2. Load + combine ----
cat("[2/6] Loading RDS...\n")
cells_list <- lapply(ckpts, function(f) {
  tryCatch(readRDS(f), error = function(e) NULL)
})
cells_list <- cells_list[!sapply(cells_list, is.null)]
cells_df <- bind_rows(cells_list)
cat("  Total cells:", nrow(cells_df), "\n")
cat("  Columns:", paste(colnames(cells_df), collapse = ", "), "\n")
stopifnot(nrow(cells_df) > 0)

# ---- 3. Filter EAE + 3 timepoints ----
cat("[3/6] Filtering EAE and Early/Peak/Late...\n")
genes <- c("B4galt6", "Ugcg", "Smpd3", "Sgpl1", "S1pr1")
stopifnot(all(genes %in% colnames(cells_df)))
stopifnot(all(c("orig.ident","timepoint","condition") %in% colnames(cells_df)))

cells_eae <- cells_df %>%
  filter(condition == "EAE",
         timepoint %in% c("Early", "Peak", "Late"))
cat("  EAE cells:", nrow(cells_eae), "\n")
cat("  Mice:", length(unique(cells_eae$orig.ident)), "\n")
cat("  Per timepoint:\n"); print(table(cells_eae$timepoint))

# ---- 4. Long + positivity per mouse ----
cat("[4/6] Computing per-mouse positivity...\n")
per_mouse <- cells_eae %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "CPM") %>%
  mutate(positive = CPM > 0) %>%
  group_by(orig.ident, timepoint, gene) %>%
  summarise(
    n_cells      = n(),
    n_positive   = sum(positive),
    pct_positive = n_positive / n_cells * 100,
    mean_pos_CPM = ifelse(n_positive > 0, mean(CPM[positive]), 0),
    .groups = "drop"
  )
cat("  per_mouse rows:", nrow(per_mouse), " (expected 5 genes x 17 mice = 85)\n")
stopifnot(nrow(per_mouse) == 85)

# ---- 5. Wilcoxon per gene: Early vs Peak ----
cat("[5/6] Wilcoxon per-mouse (Early vs Peak) for 5 genes...\n")
set.seed(42)
results_list <- list()
for (g in genes) {
  d <- per_mouse %>% filter(gene == g)
  e <- d$pct_positive[d$timepoint == "Early"]
  p <- d$pct_positive[d$timepoint == "Peak"]
  w <- tryCatch(
    wilcox.test(e, p, exact = TRUE),
    error = function(err) {
      cat("    [warn]", g, "exact failed, using normal:", conditionMessage(err), "\n")
      wilcox.test(e, p, exact = FALSE)
    }
  )
  results_list[[g]] <- data.frame(
    gene          = g,
    n_early       = length(e),
    n_peak        = length(p),
    mean_pct_early = mean(e),
    mean_pct_peak  = mean(p),
    delta_pct     = mean(p) - mean(e),
    p_raw         = w$p.value,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %-8s: Early=%.2f%% Peak=%.2f%% delta=%+.2f%% p=%.3e\n",
              g, mean(e), mean(p), mean(p) - mean(e), w$p.value))
}
results <- bind_rows(results_list)
results$p_adj_BH <- p.adjust(results$p_raw, method = "BH")
results$sig <- ifelse(results$p_adj_BH < 0.001, "***",
                ifelse(results$p_adj_BH < 0.01,  "**",
                ifelse(results$p_adj_BH < 0.05,  "*", "ns")))

# ---- 6. Plot data (mean ± SD per gene x timepoint) ----
cat("[6/6] Writing outputs...\n")
plot_data <- per_mouse %>%
  group_by(gene, timepoint) %>%
  summarise(
    mean_pct   = mean(pct_positive),
    sd_pct     = sd(pct_positive),
    n_mice     = n(),
    mean_cpm   = mean(mean_pos_CPM),
    sd_cpm     = sd(mean_pos_CPM),
    .groups = "drop"
  ) %>%
  left_join(results %>% select(gene, p_adj_BH, sig), by = "gene")

write.csv(results,    "results/snrnaseq/36c_wilcoxon_per_mouse.csv",      row.names = FALSE)
write.csv(per_mouse,  "results/snrnaseq/36c_positivity_per_mouse.csv",    row.names = FALSE)
write.csv(plot_data,  "results/snrnaseq/36c_positivity_plot_data.csv",    row.names = FALSE)

cat("\n[OK] 3 CSVs written\n\n")
cat("══════════════════════════════════════════════════\n")
cat("RESULTS (per-mouse Wilcoxon)\n")
cat("══════════════════════════════════════════════════\n")
print(results[, c("gene","n_early","n_peak","mean_pct_early","mean_pct_peak",
                  "delta_pct","p_raw","p_adj_BH","sig")], row.names = FALSE)
cat("\n══════════════════════════════════════════════════\n")
cat("PLOT DATA (mean ± SD)\n")
cat("══════════════════════════════════════════════════\n")
print(as.data.frame(plot_data), row.names = FALSE, digits = 3)

cat("\n=== FIN ===\n")
