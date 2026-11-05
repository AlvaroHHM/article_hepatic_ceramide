# ============================================================
# 36b_damage_interferon_tests.R
# 3 primary tests + 6 exploratory tests for the Damage/Interferon
# dissociation in MOL3_EAE (Figure 3C).
#
# Primary tier (caption citation):
#   1. Damage Late vs Early    (Wilcoxon exact, n=6 vs 6)
#   2. Interferon Peak vs Late (Wilcoxon exact, n=5 vs 6)
#   3. Kruskal-Wallis on delta = Damage - Interferon across stages
#
# Exploratory tier (supplementary only, not cited in caption):
#   - Damage Early vs Peak, Damage Peak vs Late
#   - Interferon Early vs Peak, Interferon Early vs Late
#   - Paired Damage vs Interferon within Late
#   - lmer interaction (marker x stage, sensitivity check)
#
# Output:
#   results/snrnaseq/36b_damage_interferon_pvals.csv
#   logs/snrnaseq/36b_damage_interferon_tests_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

dir.create("logs/snrnaseq", recursive = TRUE, showWarnings = FALSE)
sink("logs/snrnaseq/36b_damage_interferon_tests_log.txt", split = TRUE)
on.exit(sink())

cat("=== 36b: Damage/Interferon dissociation tests ===\n")
cat("Date:", format(Sys.time()), "\n\n")

ANN <- "../data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
stopifnot(file.exists(ANN))
OUT_CSV <- "results/snrnaseq/36b_damage_interferon_pvals.csv"

cat("[1/7] Loading annotation...\n")
ann <- fread(ANN, data.table = FALSE, check.names = FALSE)
stopifnot("cellTypeFalcao2018" %in% colnames(ann))
stopifnot(all(c("orig.ident","DamageScore","InterferonScore") %in% colnames(ann)))
cat("  Total cells:", nrow(ann), "\n")

cat("[2/7] Filtering MOL3_EAE...\n")
mol3 <- ann[ann$cellTypeFalcao2018 == "MOL3_EAE",
            c("orig.ident","DamageScore","InterferonScore")]
cat("  MOL3_EAE cells:", nrow(mol3), "\n")
stopifnot(nrow(mol3) > 0)

cat("[3/7] Deriving stage and filtering EAE...\n")
mol3$stage <- sub("^(Early|Peak|Late).*$", "\\1", mol3$orig.ident)
mol3$condition <- ifelse(grepl("Ctrl", mol3$orig.ident), "Ctrl", "EAE")
mol3_eae <- mol3[mol3$condition == "EAE" & mol3$stage %in% c("Early","Peak","Late"), ]
cat("  EAE cells:", nrow(mol3_eae), " | mice:", length(unique(mol3_eae$orig.ident)), "\n")
cat("  Per stage:\n"); print(table(mol3_eae$stage))
stopifnot(all(c("Early","Peak","Late") %in% mol3_eae$stage))

cat("\n[4/7] Pseudobulk per mouse...\n")
pseudo <- mol3_eae %>%
  group_by(orig.ident, stage) %>%
  summarise(
    mean_Damage    = mean(DamageScore, na.rm = TRUE),
    mean_Interferon = mean(InterferonScore, na.rm = TRUE),
    n_cells        = n(),
    .groups = "drop"
  ) %>%
  mutate(delta = mean_Damage - mean_Interferon)
cat("  Pseudobulk rows:", nrow(pseudo), " (expected 17)\n")
print(as.data.frame(pseudo))

wilcox_safe <- function(x, y, paired = FALSE) {
  tryCatch(
    wilcox.test(x, y, exact = TRUE, paired = paired),
    error = function(e) {
      cat("    [warn] exact test failed, using normal approximation:", conditionMessage(e), "\n")
      wilcox.test(x, y, exact = FALSE, paired = paired)
    }
  )
}

cat("\n[5/7] PRIMARY tests (caption citation)...\n")

d_late  <- pseudo$mean_Damage[pseudo$stage == "Late"]
d_early <- pseudo$mean_Damage[pseudo$stage == "Early"]
d_peak  <- pseudo$mean_Damage[pseudo$stage == "Peak"]
i_late  <- pseudo$mean_Interferon[pseudo$stage == "Late"]
i_early <- pseudo$mean_Interferon[pseudo$stage == "Early"]
i_peak  <- pseudo$mean_Interferon[pseudo$stage == "Peak"]

w_d_late_early <- wilcox_safe(d_late, d_early)
cat(sprintf("  P1 Damage_Late_vs_Early:    n=%d vs %d | mean=%.4f vs %.4f | p=%.4f\n",
            length(d_late), length(d_early), mean(d_late), mean(d_early), w_d_late_early$p.value))

w_i_peak_late <- wilcox_safe(i_peak, i_late)
cat(sprintf("  P2 Interferon_Peak_vs_Late: n=%d vs %d | mean=%.4f vs %.4f | p=%.4f\n",
            length(i_peak), length(i_late), mean(i_peak), mean(i_late), w_i_peak_late$p.value))

kw_delta <- kruskal.test(delta ~ stage, data = pseudo)
cat(sprintf("  P3 KW_delta (across stages): H=%.3f, df=%d, p=%.4f\n",
            kw_delta$statistic, kw_delta$parameter, kw_delta$p.value))

p_primary <- c(w_d_late_early$p.value, w_i_peak_late$p.value)
p_primary_adj <- p.adjust(p_primary, method = "BH")

cat("\n[6/7] EXPLORATORY tests (supplementary only)...\n")

w_d_e_p <- wilcox_safe(d_early, d_peak)
w_d_p_l <- wilcox_safe(d_peak, d_late)
w_i_e_p <- wilcox_safe(i_early, i_peak)
w_i_e_l <- wilcox_safe(i_early, i_late)
w_paired_late <- wilcox_safe(d_late, i_late, paired = TRUE)

cat(sprintf("  E1 Damage_Early_vs_Peak:             p=%.4f\n", w_d_e_p$p.value))
cat(sprintf("  E2 Damage_Peak_vs_Late:              p=%.4f\n", w_d_p_l$p.value))
cat(sprintf("  E3 Interferon_Early_vs_Peak:         p=%.4f\n", w_i_e_p$p.value))
cat(sprintf("  E4 Interferon_Early_vs_Late:         p=%.4f\n", w_i_e_l$p.value))
cat(sprintf("  E5 Paired_Damage_vs_Interferon_Late: p=%.4f\n", w_paired_late$p.value))

p_explor <- c(w_d_e_p$p.value, w_d_p_l$p.value, w_i_e_p$p.value, w_i_e_l$p.value, w_paired_late$p.value)
p_explor_adj <- p.adjust(p_explor, method = "BH")

long <- rbind(
  data.frame(mouse = pseudo$orig.ident, stage = pseudo$stage,
             marker = "Damage",     score = pseudo$mean_Damage),
  data.frame(mouse = pseudo$orig.ident, stage = pseudo$stage,
             marker = "Interferon", score = pseudo$mean_Interferon)
)
long$stage  <- factor(long$stage,  levels = c("Early","Peak","Late"))
long$marker <- factor(long$marker, levels = c("Damage","Interferon"))

lmer_result <- tryCatch({
  suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
  m <- lmer(score ~ marker * stage + (1 | mouse), data = long)
  aov_tab <- anova(m)
  list(p    = aov_tab["marker:stage", "Pr(>F)"],
       F    = aov_tab["marker:stage", "F value"],
       note = "lmer_converged")
}, error = function(e) {
  list(p = NA, F = NA, note = paste0("lmer_failed: ", conditionMessage(e)))
})
cat(sprintf("  E6 lmer_marker_x_stage: F=%s, p=%s (%s)\n",
            ifelse(is.na(lmer_result$F), "NA", sprintf("%.3f", lmer_result$F)),
            ifelse(is.na(lmer_result$p), "NA", sprintf("%.4f", lmer_result$p)),
            lmer_result$note))

cat("\n[7/7] Writing CSV...\n")

results <- data.frame(
  test    = c("Damage_Late_vs_Early","Interferon_Peak_vs_Late","KW_delta",
              "Damage_Early_vs_Peak","Damage_Peak_vs_Late",
              "Interferon_Early_vs_Peak","Interferon_Early_vs_Late",
              "Paired_Damage_vs_Interferon_Late","lmer_marker_x_stage"),
  tier    = c("primary","primary","primary",
              rep("exploratory", 6)),
  marker  = c("Damage","Interferon","combined",
              "Damage","Damage","Interferon","Interferon","combined","combined"),
  comparison = c("Late_vs_Early","Peak_vs_Late","delta_across_stages",
                 "Early_vs_Peak","Peak_vs_Late",
                 "Early_vs_Peak","Early_vs_Late","paired_Late","marker_x_stage"),
  n1 = c(length(d_late), length(i_peak), nrow(pseudo),
         length(d_early), length(d_peak),
         length(i_early), length(i_early),
         length(d_late), nrow(long)),
  n2 = c(length(d_early), length(i_late), NA,
         length(d_peak), length(d_late),
         length(i_peak), length(i_late),
         length(i_late), NA),
  mean1 = c(mean(d_late), mean(i_peak), NA,
            mean(d_early), mean(d_peak),
            mean(i_early), mean(i_early),
            mean(d_late), NA),
  mean2 = c(mean(d_early), mean(i_late), NA,
            mean(d_peak), mean(d_late),
            mean(i_peak), mean(i_late),
            mean(i_late), NA),
  p_raw = c(w_d_late_early$p.value, w_i_peak_late$p.value, kw_delta$p.value,
            w_d_e_p$p.value, w_d_p_l$p.value, w_i_e_p$p.value, w_i_e_l$p.value,
            w_paired_late$p.value, lmer_result$p),
  p_adj = c(p_primary_adj, NA, p_explor_adj, NA),
  note  = c(rep("exact_Wilcoxon", 2), "KW_H_2",
            rep("exact_Wilcoxon", 5),
            lmer_result$note)
)

sig_star <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05, "*", "ns"))))
}
results$sig <- sig_star(results$p_adj)

write.csv(results, OUT_CSV, row.names = FALSE)
cat("[OK]", OUT_CSV, "\n\n")

cat("===============================================================\n")
cat("PRIMARY (cited in caption)\n")
cat("===============================================================\n")
print(results[results$tier == "primary",
              c("test","n1","n2","mean1","mean2","p_raw","p_adj","sig","note")],
      row.names = FALSE)

cat("\n===============================================================\n")
cat("EXPLORATORY (supplementary only)\n")
cat("===============================================================\n")
print(results[results$tier == "exploratory",
              c("test","n1","n2","mean1","mean2","p_raw","p_adj","sig","note")],
      row.names = FALSE)

cat("\n===============================================================\n")
cat("DONE\n")
cat("===============================================================\n")
