# ============================================================
# 36d_glmer_per_mouse.R
# glmer binomial per-mouse (with random effect) + t-test logit fallback
# for the 5 Figure 5 genes. Compares with Wilcoxon per-mouse (36c).
#
# Motivation: Wilcoxon per-mouse uses only % per mouse (n=6 vs 5),
# losing the n_cells information. glmer binomial uses cbind(n_pos, n_neg)
# per mouse + random effect, potentially with more power.
#
# Outputs:
#   results/snrnaseq/36d_glmer_per_mouse.csv
#   logs/snrnaseq/36d_glmer_per_mouse_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(lmerTest)
})
dir.create("logs/snrnaseq", recursive = TRUE, showWarnings = FALSE)
sink("logs/snrnaseq/36d_glmer_per_mouse_log.txt", split = TRUE)
on.exit(sink())

cat("=== 36d: glmer binomial per-mouse (5 genes) ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# ---- 1. Checkpoints ----
ckpt_candidates <- c("results/35b_checkpoint", "../results/35b_checkpoint")
ckpt_dir <- ckpt_candidates[dir.exists(ckpt_candidates)][1]
if (is.na(ckpt_dir)) stop("No checkpoint dir found")
ckpts <- list.files(ckpt_dir, pattern = "\\.rds$", full.names = TRUE)
cat("[1/5] Checkpoints:", length(ckpts), "\n")

# ---- 2. Load + filter EAE ----
cells_list <- lapply(ckpts, function(f) tryCatch(readRDS(f), error = function(e) NULL))
cells_list <- cells_list[!sapply(cells_list, is.null)]
cells_df <- bind_rows(cells_list)
cat("[2/5] Total cells:", nrow(cells_df), "\n")

genes <- c("B4galt6", "Ugcg", "Smpd3", "Sgpl1", "S1pr1")
cells_eae <- cells_df %>%
  filter(condition == "EAE", timepoint %in% c("Early", "Peak", "Late"))
cat("  EAE cells:", nrow(cells_eae), " | mice:", length(unique(cells_eae$orig.ident)), "\n")

# ---- 3. Per-mouse n_pos / n_neg per gene ----
cat("[3/5] Per-mouse binomial data...\n")
per_mouse <- cells_eae %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "CPM") %>%
  mutate(positive = CPM > 0) %>%
  group_by(orig.ident, timepoint, gene) %>%
  summarise(n_pos = sum(positive), n_neg = sum(!positive), .groups = "drop") %>%
  mutate(stage = factor(timepoint, levels = c("Early", "Peak", "Late")))

# ---- 4. glmer per gene + logit t-test fallback ----
cat("[4/5] Fitting models per gene...\n")
results_list <- list()
for (g in genes) {
  d <- per_mouse %>% filter(gene == g)
  # A) glmer binomial with random effect (may not converge with n=6v5)
  glmer_res <- tryCatch({
    m <- glmer(cbind(n_pos, n_neg) ~ stage + (1 | orig.ident),
               data = d, family = binomial,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    coefs <- summary(m)$coefficients
    if ("stagePeak" %in% rownames(coefs)) {
      list(p = coefs["stagePeak", "Pr(>|z|)"],
           coef = coefs["stagePeak", "Estimate"],
           note = "glmer_converged")
    } else list(p = NA, coef = NA, note = "glmer_no_stagePeak")
  }, error = function(e) list(p = NA, coef = NA, note = paste0("glmer_failed: ", conditionMessage(e))))

  # B) t-test on logit-transformed proportions (stable fallback)
  d_wide <- d %>%
    mutate(pct = n_pos / (n_pos + n_neg),
           logit = log((pct + 0.001) / (1 - pct + 0.001))) %>%
    select(orig.ident, stage, logit) %>%
    pivot_wider(names_from = stage, values_from = logit)
  t_res <- tryCatch({
    e <- d_wide$Early; pk <- d_wide$Peak
    tt <- t.test(pk, e)
    list(p = tt$p.value, coef = mean(pk) - mean(e), note = "t_test_logit")
  }, error = function(e) list(p = NA, coef = NA, note = paste0("t_failed: ", conditionMessage(e))))

  results_list[[g]] <- data.frame(
    gene = g,
    n_early = sum(d$timepoint == "Early"),
    n_peak  = sum(d$timepoint == "Peak"),
    glmer_p = glmer_res$p, glmer_coef = glmer_res$coef, glmer_note = glmer_res$note,
    ttest_logit_p = t_res$p, ttest_logit_coef = t_res$coef, ttest_note = t_res$note,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %-8s: glmer p=%s (%s) | t-logit p=%.4f\n",
              g,
              ifelse(is.na(glmer_res$p), "NA", sprintf("%.4e", glmer_res$p)),
              glmer_res$note, t_res$p))
}
results <- bind_rows(results_list)
results$glmer_p_adj <- p.adjust(results$glmer_p, "BH")
results$ttest_logit_p_adj <- p.adjust(results$ttest_logit_p, "BH")

# ---- 5. Write ----
cat("[5/5] Writing output...\n")
write.csv(results, "results/snrnaseq/36d_glmer_per_mouse.csv", row.names = FALSE)

cat("\n══════════════════════════════════════════════════════\n")
cat("RESULTS — glmer binomial vs t-test logit (per-mouse)\n")
cat("══════════════════════════════════════════════════════\n")
print(results[, c("gene","n_early","n_peak","glmer_p","glmer_p_adj","glmer_note",
                  "ttest_logit_p","ttest_logit_p_adj")], row.names = FALSE)
cat("\n=== FIN ===\n")
