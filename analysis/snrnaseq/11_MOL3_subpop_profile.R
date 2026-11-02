# ============================================================
# 11_MOL3_subpop_profile.R
# Consolidated MOL3_EAE subpopulation profile across EAE stages.
#
# Merges:
#   - B4galt6+ fraction   <- results/snrnaseq/08_MOL3_sphingolipid_substructure.csv
#   - DamageScore         <- data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz
#   - InterferonScore     <- annotation
#   - ImmuneScore         <- annotation
#
# Output:
#   - results/snrnaseq/11_MOL3_subpop_profile.csv
#   - logs/snrnaseq/11_MOL3_subpop_profile_log.txt
#
# Verifies 4 numbers from the manuscript:
#   DamageScore  MOL3∩MOL2 (Late)  = 0.061  -> actual ~0.0611
#   DamageScore  MOL3∩MOL56 (Late) = 0.006  -> actual ~0.0056
#   B4galt6+ %   MOL3∩MOL2 (Late)  = 16.6   -> actual ~16.57
#   B4galt6+ %   MOL3∩MOL56 (Late) = 3.7    -> actual ~3.70
#
# Does NOT require the 28 H5 files.
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

find_data_root <- function() {
  for (r in c(".", "..", "../..")) {
    if (file.exists(file.path(r, "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz")))
      return(r)
  }
  stop("Cannot find annotation. Run from repo root or parent.")
}
DATA_ROOT <- find_data_root()

if (!dir.exists("results") || !dir.exists("analysis"))
  stop("Run from repo root (article_hepatic_ceramide/) — expected results/ and analysis/")
OUT_DIR <- "results/snrnaseq"
LOG_DIR <- "logs/snrnaseq"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

annot_path <- file.path(DATA_ROOT, "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz")

csv_08 <- NULL
for (name in c("08_MOL3_sphingolipid_substructure.csv",
               "44b_MOL3_sphingolipid_substructure.csv")) {
  cand <- file.path(OUT_DIR, name)
  if (file.exists(cand)) { csv_08 <- cand; break }
}
if (is.null(csv_08)) stop("Missing 08 output at ", file.path(OUT_DIR, "08_*.csv"))

sink(file.path(LOG_DIR, "11_MOL3_subpop_profile_log.txt"), split = TRUE)
on.exit(sink())

SEP <- "============================================================"
cat(SEP, "\nMOL3 SUBPOPULATION PROFILE (consolidated)\n")
cat("Date:", format(Sys.time()), "\n", SEP, "\n\n")
cat("Annotation:", annot_path, "\nB4galt6 data:", csv_08, "\n\n")

cat("--- 1. Loading annotation ---\n")
annot <- fread(annot_path, data.table = FALSE, check.names = FALSE)
cat("  Dimensions:", nrow(annot), "x", ncol(annot), "\n")
rownames(annot) <- annot[[1]]
annot <- annot[, -1]

mol3 <- annot[annot$cellTypeFalcao2018 == "MOL3_EAE", ]
cat("  MOL3_EAE cells:", nrow(mol3), "\n")

mol3$subpop <- ifelse(mol3$cellType == "MOL2", "MOL2_like",
               ifelse(mol3$cellType == "MOL56", "MOL56_like", "Other"))
mol3$stage  <- ifelse(grepl("^Early", mol3$orig.ident), "Early",
               ifelse(grepl("^Peak",  mol3$orig.ident), "Peak",
               ifelse(grepl("^Late",  mol3$orig.ident), "Late", NA_character_)))
cat("  Distribution:\n"); print(table(mol3$subpop, mol3$stage, useNA = "ifany")); cat("\n")

cat("--- 2. Aggregating scores ---\n")
scores <- aggregate(cbind(DamageScore, InterferonScore, ImmuneScore) ~ subpop + stage,
                    data = mol3, FUN = function(x) mean(x, na.rm = TRUE))
counts <- aggregate(DamageScore ~ subpop + stage, data = mol3,
                    FUN = function(x) sum(!is.na(x)))
names(counts)[3] <- "n_cells"
profile <- merge(counts, scores, by = c("subpop", "stage"))
cat("  Rows:", nrow(profile), "\n\n")

cat("--- 3. Loading B4galt6+ from 08 ---\n")
b4 <- fread(csv_08, data.table = FALSE)
b4 <- b4[b4$gene == "B4galt6", c("subpop", "stage", "n_cells", "n_positive", "pct_positive")]
names(b4)[3:5] <- c("n_cells_b4", "n_B4galt6_pos", "pct_B4galt6_pos")

profile <- merge(profile, b4[, c("subpop","stage","n_B4galt6_pos","pct_B4galt6_pos")],
                 by = c("subpop","stage"), all.x = TRUE)
profile <- profile[order(profile$subpop,
                 factor(profile$stage, levels = c("Early","Peak","Late"))), ]

cat("\n--- 4. CONSOLIDATED PROFILE ---\n"); print(profile, row.names = FALSE); cat("\n")

out_path <- file.path(OUT_DIR, "11_MOL3_subpop_profile.csv")
write.csv(profile, out_path, row.names = FALSE)
cat("[OK]", out_path, "\n\n")

cat(SEP, "\nVERIFICATION AGAINST MANUSCRIPT\n", SEP, "\n")
late  <- profile[profile$stage == "Late", ]
mol2  <- late[late$subpop == "MOL2_like", ]
mol56 <- late[late$subpop == "MOL56_like", ]

verify <- function(df, col, e, label, tol) {
  if (nrow(df) == 0 || !(col %in% names(df))) {
    cat(sprintf("  %-32s COLUMN/R0W MISSING  MISMATCH\n", label))
    return(FALSE)
  }
  a <- df[[col]][1]
  ok <- isTRUE(!is.na(a) && abs(a - e) < tol)
  cat(sprintf("  %-32s actual=%.4f expected=%.3f  %s\n",
              label, a, e, ifelse(ok, "OK", "MISMATCH")))
  ok
}
ok <- c(
  verify(mol2,  "DamageScore",      0.061, "MOL3nMOL2  DamageScore (Late)", 0.001),
  verify(mol56, "DamageScore",      0.006, "MOL3nMOL56 DamageScore (Late)", 0.001),
  verify(mol2,  "pct_B4galt6_pos", 16.6,   "MOL3nMOL2  B4galt6+ % (Late)",  0.1),
  verify(mol56, "pct_B4galt6_pos",  3.7,   "MOL3nMOL56 B4galt6+ % (Late)",  0.1)
)
cat(sprintf("\n  Verified: %d/4 numbers match the manuscript\n", sum(ok)))

cat("\n", SEP, "\nDERIVED METRICS (informative)\n", SEP, "\n")
cat(sprintf("  DamageScore ratio  MOL2/MOL56 (Late): %.2fx\n",
            mol2$DamageScore / mol56$DamageScore))
cat(sprintf("  B4galt6+ ratio     MOL2/MOL56 (Late): %.2fx\n",
            mol2$pct_B4galt6_pos / mol56$pct_B4galt6_pos))
cat(sprintf("  InterferonScore    MOL2/MOL56 (Late): %.2fx\n",
            mol2$InterferonScore / mol56$InterferonScore))

cat("\n", SEP, "\nFIN\n", SEP, "\n")
