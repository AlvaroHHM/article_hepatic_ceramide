# ============================================================
# 31_verify_v8_numbers.R
# Verifica los numeros citados en ANTEPROYECTO_HOLOBIONTE_v8.tex
# contra los CSVs fuente. Genera PASS/WARN/FAIL.
#
# Uso: Rscript scripts/31_verify_v8_numbers.R
# Salidas:
#   results/31_verification_results.csv
#   results/31_n_per_timepoint.csv  (tabla que falta en v8)
#   logs/31_verify_v8_numbers_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

sink("logs/31_verify_v8_numbers_log.txt", split = TRUE)
on.exit(sink())

cat("============================================================\n")
cat("VERIFICATION REPORT: ANTEPROYECTO_HOLOBIONTE_v8\n")
cat("Date:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# --- Helper: flexible file finder ---
find_file <- function(candidates) {
  for (f in candidates) if (file.exists(f)) return(f)
  return(NA_character_)
}

# --- Helper: PASS/WARN/FAIL table ---
results_tbl <- data.frame(
  check    = character(),
  item     = character(),
  src_val  = numeric(),
  paper_val = numeric(),
  delta    = numeric(),
  status   = character(),
  stringsAsFactors = FALSE
)

add_check <- function(check, item, src, paper, tol = 0.05) {
  d <- ifelse(is.na(src) | is.na(paper), NA_real_, src - paper)
  st <- ifelse(is.na(d), "MISSING",
        ifelse(abs(d) <= tol, "PASS",
        ifelse(abs(d) <= 3 * tol, "WARN", "FAIL")))
  results_tbl[nrow(results_tbl) + 1, ] <<- list(
    check = check, item = item,
    src_val = src, paper_val = paper,
    delta = d, status = st
  )
  cat(sprintf("  %-42s  src=%-10.4f  paper=%-10.4f  d=%+.4f  [%s]\n",
              item, src, paper, d, st))
}

# --- Locate source files ---
falcao_file   <- find_file(c("results/30b_falcao_doble_jab_eae_only.csv",
                             "results/30b_falcao_doble_jab_table.csv"))
mapping_file  <- find_file(c("results/34_GSM_mapping_v2.csv",
                             "results/34_GSM_mapping.csv"))
gene_file     <- find_file(c("results/35b_MOL3_early_vs_peak.csv",
                             "results/35b_MOL3_early_vs_peak_clean.csv"))
damage_file   <- find_file(c("results/36_damage_mol3_by_timepoint.csv"))
positivity_f  <- find_file(c("results/36_positivity.csv"))

cat("Sources:\n")
cat("  falcao:    ", falcao_file,  "\n", sep = "")
cat("  mapping:   ", mapping_file, "\n", sep = "")
cat("  gene:      ", gene_file,    "\n", sep = "")
cat("  damage:    ", damage_file,  "\n", sep = "")
cat("  positivity:", positivity_f, "\n\n", sep = "")

# ============================================================
# CHECK 1: Falcao VI and PI
# ============================================================
cat("--- CHECK 1: Falcao VI/PI ---\n")
if (!is.na(falcao_file)) {
  falcao <- read.csv(falcao_file, check.names = FALSE)
  cat("  Columns:", paste(colnames(falcao), collapse = ", "), "\n")
  for (i in seq_len(nrow(falcao))) {
    st <- falcao$state[i]
    e  <- falcao$Early_EAE[i]
    p  <- falcao$Peak_EAE[i]
    l  <- falcao$Late_EAE[i]
    vi_calc <- if (e > 0) (p - e) / e else NA_real_
    pi_calc <- if (p > 0) l / p else NA_real_
    add_check("1_VI", paste0(st, " VI"), vi_calc, falcao$VI[i], tol = 0.05)
    add_check("1_PI", paste0(st, " PI"), pi_calc, falcao$PI[i], tol = 0.02)
  }
} else cat("  SKIP: falcao file missing\n")

# ============================================================
# CHECK 2: tab:mol3-share internal consistency
# ============================================================
cat("\n--- CHECK 2: tab:mol3-share (ratios) ---\n")
share <- data.frame(
  tp       = c("Early", "Peak", "Late"),
  mol3     = c(4.9,     20.6,    28.9),
  eae_tot  = c(65.0,    92.3,    48.5),
  ratio_p  = c(7.5,     22.4,    59.6)
)
for (i in seq_len(nrow(share))) {
  r <- share$mol3[i] / share$eae_tot[i] * 100
  add_check("2_share", paste0(share$tp[i], " mol3/eae ratio"),
            r, share$ratio_p[i], tol = 0.15)
}

# ============================================================
# CHECK 3: text claims (71.7, 48.5, 92.3)
# ============================================================
cat("\n--- CHECK 3: text claims ---\n")
add_check("3_text", "71.7 = 92.3 - 20.6", 92.3 - 20.6, 71.7, tol = 0.05)
add_check("3_text", "48.5 = 28.9/0.596",  28.9 / 0.596, 48.5, tol = 0.10)
add_check("3_text", "92.3 = 20.6/0.224",  20.6 / 0.224, 92.3, tol = 0.30)

# ============================================================
# CHECK 4: EAE-only totals
# ============================================================
cat("\n--- CHECK 4: EAE-only totals (sum of states) ---\n")
if (!is.na(falcao_file)) {
  e_sum <- sum(falcao$Early_EAE) * 100
  p_sum <- sum(falcao$Peak_EAE)  * 100
  l_sum <- sum(falcao$Late_EAE)  * 100
  add_check("4_sum", "Early EAE-OLG total", e_sum, 65.0, tol = 0.5)
  add_check("4_sum", "Peak EAE-OLG total",  p_sum, 92.3, tol = 0.3)
  add_check("4_sum", "Late EAE-OLG total",  l_sum, 48.5, tol = 0.3)
} else cat("  SKIP: falcao file missing\n")

# ============================================================
# CHECK 5: Late composition sums to 100
# ============================================================
cat("\n--- CHECK 5: Late composition sums to 100 ---\n")
late_parts <- c(EAE_OLG = 48.5, Homeostatic = 43.6, Progenitor = 7.8,
                Microglia = 0.1)
add_check("5_sum100", "Late parts sum", sum(late_parts), 100.0, tol = 0.2)

# ============================================================
# CHECK 6: n per timepoint (mapping file) + write table
# ============================================================
cat("\n--- CHECK 6: n per timepoint ---\n")
if (!is.na(mapping_file)) {
  mp <- read.csv(mapping_file, check.names = FALSE)
  cat("  Columns:", paste(colnames(mp), collapse = ", "), "\n")
  # Try to detect columns
  col_cond <- intersect(c("condition","Condition","group","Group"),
                        colnames(mp))[1]
  col_time <- intersect(c("timepoint","time_point","tp","Timepoint",
                          "time","Time"),
                        colnames(mp))[1]
  if (!is.na(col_cond) && !is.na(col_time)) {
    tab <- mp %>%
      count(.data[[col_cond]], .data[[col_time]]) %>%
      tidyr::pivot_wider(names_from = .data[[col_cond]],
                         values_from = n, values_fill = 0)
    cat("\n  Sample counts (condition x timepoint):\n")
    print(as.data.frame(tab))
    write.csv(tab, "results/31_n_per_timepoint.csv", row.names = FALSE)
    cat("  -> results/31_n_per_timepoint.csv\n")
    # Cross-check paper claims
    for (tp in c("Early","Peak","Late")) {
      tp_upper <- toupper(tp)
      row <- tab[grepl(tp_upper, toupper(as.character(tab[[col_time]])),
                       ignore.case = TRUE), ]
      if (nrow(row) > 0) {
        eae_n <- sum(row[, grep("EAE", colnames(row), value = TRUE),
                         drop = TRUE], na.rm = TRUE)
        cat(sprintf("  %-6s EAE samples: %d\n", tp, eae_n))
      }
    }
    # Total counts
    tot <- table(mp[[col_cond]])
    cat("\n  Total samples per condition:\n")
    print(tot)
  } else {
    cat("  SKIP: could not detect condition/time columns\n")
  }
} else {
  cat("  SKIP: mapping file missing\n")
  cat("  Paper claims: Early=6 EAE, Peak=5 EAE, Late=6 EAE, 17 EAE + 9 Ctrl = 26\n")
}

# ============================================================
# CHECK 7: gene expression (35b)
# ============================================================
cat("\n--- CHECK 7: gene expression Early->Peak ---\n")
if (!is.na(gene_file)) {
  ge <- read.csv(gene_file, check.names = FALSE)
  cat("  Columns:", paste(colnames(ge), collapse = ", "), "\n")
  # Paper targets
  paper <- list(
    B4galt6 = c(3.78, 10.1, +1.22),
    Ugcg    = c(4.43,  8.68, +0.83),
    Smpd3   = c(2.92,  6.05, +0.85),
    Sgpl1   = c(5.05,  9.90, +0.85),
    Sgms1   = c(25.7,  39.5, +0.60),
    S1pr1   = c(1.64,  0.25, -1.08)
  )
  for (g in names(paper)) {
    row <- ge[ge$gene == g, ]
    if (nrow(row) == 0) {
      cat(sprintf("  %-8s: NOT IN FILE\n", g))
      next
    }
    # Try to find the relevant columns
    e_col <- grep("early", colnames(ge), ignore.case = TRUE, value = TRUE)[1]
    p_col <- grep("peak",  colnames(ge), ignore.case = TRUE, value = TRUE)[1]
    fc_col <- grep("log2FC", colnames(ge), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(e_col) || is.na(p_col) || is.na(fc_col)) {
      cat(sprintf("  %-8s: columns not detected\n", g)); next
    }
    add_check("7_gene", paste0(g, " Early"), row[[e_col]], paper[[g]][1], tol = 0.5)
    add_check("7_gene", paste0(g, " Peak"),  row[[p_col]], paper[[g]][2], tol = 0.5)
    add_check("7_gene", paste0(g, " log2FC"), row[[fc_col]], paper[[g]][3], tol = 0.1)
  }
} else cat("  SKIP: gene file missing\n")

# ============================================================
# CHECK 8: damage/immune/interferon scores
# ============================================================
cat("\n--- CHECK 8: damage scores ---\n")
if (!is.na(damage_file)) {
  dm <- read.csv(damage_file, check.names = FALSE)
  cat("  Columns:", paste(colnames(dm), collapse = ", "), "\n")
  print(dm)
  # Paper: Early D=-0.002 I=0.064 F=0.343
  #        Peak  D=+0.072 I=0.067 F=0.562
  #        Late  D=+0.038 I=0.061 F=0.112
} else cat("  SKIP: damage file missing\n")

# ============================================================
# SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")
if (nrow(results_tbl) > 0) {
  cnt <- table(factor(results_tbl$status,
                     levels = c("PASS","WARN","FAIL","MISSING")))
  cat(sprintf("  PASS:    %d\n", cnt["PASS"]))
  cat(sprintf("  WARN:    %d\n", cnt["WARN"]))
  cat(sprintf("  FAIL:    %d\n", cnt["FAIL"]))
  cat(sprintf("  MISSING: %d\n", cnt["MISSING"]))

  fails <- results_tbl[results_tbl$status == "FAIL", ]
  if (nrow(fails) > 0) {
    cat("\nFAILURES (require .tex edit):\n")
    for (i in seq_len(nrow(fails))) {
      cat(sprintf("  [%s] %-40s src=%.4f paper=%.4f\n",
                  fails$check[i], fails$item[i],
                  fails$src_val[i], fails$paper_val[i]))
    }
  }

  warns <- results_tbl[results_tbl$status == "WARN", ]
  if (nrow(warns) > 0) {
    cat("\nWARNINGS (acceptable rounding):\n")
    for (i in seq_len(nrow(warns))) {
      cat(sprintf("  %-40s d=%+.4f\n", warns$item[i], warns$delta[i]))
    }
  }

  write.csv(results_tbl, "results/31_verification_results.csv",
            row.names = FALSE)
  cat("\nFull results: results/31_verification_results.csv\n")
} else {
  cat("  No checks executed (all sources missing).\n")
}

cat("\n=== FIN ===\n")
