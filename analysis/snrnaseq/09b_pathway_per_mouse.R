# ============================================================
# 09b_pathway_per_mouse.R
# Per-mouse pathway module analysis for the 7 modules of
# 09_pathway_MOL2_vs_MOL56.R (replaces per-cell pseudoreplication
# with per-mouse lmer + paired t-test).
#
# Reads the 28 H5 files (~30 min) to extract all ~40 module genes
# for MOL3_EAE cells. Aggregate per mouse × subpop.
#
# Inputs:
#   data/GSE250589/h5/*.h5 (28 files, 2.3 GB)
#   data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz
#   results/snrnaseq/00_GSM_mapping.csv
# Outputs:
#   results/snrnaseq/09b_pathway_per_mouse_pvals.csv
#   results/snrnaseq/09b_per_mouse_module_scores.csv
#   results/snrnaseq/09b_percell_vs_permouse.csv
#   logs/snrnaseq/09b_pathway_per_mouse_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(Matrix);     library(hdf5r)
  library(lme4);       library(lmerTest)
})

dir.create("logs/snrnaseq", recursive = TRUE, showWarnings = FALSE)
sink("logs/snrnaseq/09b_pathway_per_mouse_log.txt", split = TRUE)
on.exit(sink())

cat("=== 09b: per-mouse pathway modules (7 modules) ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# ---- Paths ----
h5_cand <- c("data/GSE250589/h5", "../data/GSE250589/h5")
h5_dir  <- h5_cand[dir.exists(h5_cand)][1]
stopifnot(!is.na(h5_dir))
cat("[h5] ", h5_dir, "\n")

ann_cand <- c("data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz",
              "../data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz")
ann_file <- ann_cand[file.exists(ann_cand)][1]
stopifnot(!is.na(ann_file))
cat("[ann] ", ann_file, "\n")

gsm_map <- fread("results/snrnaseq/00_GSM_mapping.csv", data.table = FALSE)

# ---- Modules ----
pathways <- list(
  Sphingolipid_synthesis        = c("Sptlc1","Sptlc2","Sptlc3","Kdsr","Degs1","Degs2",
                                    "Cers1","Cers2","Cers3","Cers4","Cers5","Cers6"),
  Sphingolipid_interconversion  = c("Smpd1","Smpd2","Smpd3","Smpd4",
                                    "Sgms1","Sgms2","Asah1","Asah2",
                                    "Acer1","Acer2","Acer3","Sphk1","Sphk2","Sgpl1"),
  LacCer_synthesis              = c("Ugcg","B4galt5","B4galt6"),
  S1P_receptors                 = c("S1pr1","S1pr2","S1pr3","S1pr4","S1pr5"),
  Myelin_structural             = c("Mbp","Plp1","Mag","Mog","Cnp"),
  Apoptosis                     = c("Bax","Bak1","Casp3","Casp8","Casp9","Bcl2","Bcl2l1","Bid"),
  B4GALT6_alone                 = c("B4galt6")
)
all_genes <- unique(unlist(pathways))
cat("[modules] ", length(pathways), " modules, ", length(all_genes), " unique genes\n\n")

# ---- Annotation ----
cat("[1/9] Loading annotation...\n")
ann <- fread(ann_file, data.table = FALSE, check.names = FALSE)
rownames(ann) <- ann[[1]]; ann <- ann[, -1]
cat("  Total cells:", nrow(ann), "\n")

mol3 <- ann[ann$cellTypeFalcao2018 == "MOL3_EAE", c("orig.ident","cellType")]
mol3$barcode_clean <- sub("-.*$", "", rownames(mol3))
mol3$subpop <- ifelse(mol3$cellType == "MOL2", "MOL2_like",
               ifelse(mol3$cellType == "MOL56", "MOL56_like", "Other"))
mol3 <- mol3[mol3$subpop %in% c("MOL2_like","MOL56_like"), ]
mol3$stage <- ifelse(grepl("^Early", mol3$orig.ident), "Early",
              ifelse(grepl("^Peak",  mol3$orig.ident), "Peak",
              ifelse(grepl("^Late",  mol3$orig.ident), "Late", NA)))
mol3$condition <- ifelse(grepl("Ctrl", mol3$orig.ident), "Ctrl", "EAE")
mol3 <- mol3[mol3$condition == "EAE" & !is.na(mol3$stage), ]
mol3$GSM <- gsm_map$GSM[match(mol3$orig.ident, gsm_map$orig.ident)]
cat("  MOL3_EAE (EAE, valid stage): ", nrow(mol3), " cells | ",
    length(unique(mol3$orig.ident)), " mice\n")
cat("  Cells with GSM: ", sum(!is.na(mol3$GSM)), "\n")
print(table(mol3$subpop))

# ---- H5 reading ----
cat("\n[2/9] Reading 28 H5 files (~5-10 min)...\n")
find_h5 <- function(gsm) {
  f <- list.files(h5_dir, pattern = gsm, full.names = TRUE)
  if (length(f) == 0) NA_character_ else f[1]
}
gsms <- unique(mol3$GSM[!is.na(mol3$GSM)])
cat("  GSMs to process:", length(gsms), "\n")

read_h5_for_mol3 <- function(h5_path, gsm, target_bcs, gene_set) {
  h <- H5File$new(h5_path, mode = "r")
  on.exit(tryCatch(h$close_all(), error = function(e) NULL))
  bcs <- as.character(h[["matrix/barcodes"]][])
  g   <- as.character(h[["matrix/features/name"]][])
  sh  <- as.integer(h[["matrix/shape"]][])
  d   <- as.numeric(h[["matrix/data"]][])
  i   <- as.integer(h[["matrix/indices"]][])
  p   <- as.integer(h[["matrix/indptr"]][])
  m <- sparseMatrix(i = i + 1L, p = p, x = d, dims = sh, index1 = TRUE, giveCsparse = TRUE)
  rownames(m) <- g; colnames(m) <- bcs
  total_umi_full <- Matrix::colSums(m)

  bcs_clean    <- sub("-.*$", "", bcs)
  target_clean <- sub("-.*$", "", target_bcs)
  mi    <- match(target_clean, bcs_clean)
  valid <- !is.na(mi)
  if (sum(valid) == 0) return(list(df = NULL, n_matched = 0, n_genes = 0))
  gp <- intersect(gene_set, g)
  if (length(gp) == 0) return(list(df = NULL, n_matched = sum(valid), n_genes = 0))

  sub_counts <- m[gp, mi[valid], drop = FALSE]
  sub_umi    <- total_umi_full[mi[valid]]
  df <- as.data.frame(as.matrix(t(sub_counts)))
  colnames(df) <- gp
  df$barcode_clean <- target_bcs[valid]
  df$total_umi     <- as.numeric(sub_umi)
  df$GSM           <- gsm
  list(df = df, n_matched = sum(valid), n_genes = length(gp))
}

all_dfs <- list(); n_done <- 0
for (gsm in gsms) {
  n_done <- n_done + 1
  h5_path <- find_h5(gsm)
  if (is.na(h5_path)) {
    cat(sprintf("  [%2d/%2d] %s -> NO H5\n", n_done, length(gsms), gsm)); next
  }
  targets <- mol3[mol3$GSM == gsm, ]
  res <- tryCatch(
    read_h5_for_mol3(h5_path, gsm, targets$barcode_clean, all_genes),
    error = function(e) {
      cat(sprintf("  [%2d/%2d] %s -> ERROR: %s\n", n_done, length(gsms), gsm, conditionMessage(e))); NULL
    }
  )
  if (is.null(res) || is.null(res$df)) {
    n_m <- if (is.null(res) || is.null(res$n_matched)) 0 else res$n_matched
    n_g <- if (is.null(res) || is.null(res$n_genes))   0 else res$n_genes
    cat(sprintf("  [%2d/%2d] %s -> %d cells, %d genes (SKIP)\n",
                n_done, length(gsms), gsm, n_m, n_g)); next
  }
  cat(sprintf("  [%2d/%2d] %s -> %d cells, %d/%d genes\n",
              n_done, length(gsms), gsm, res$n_matched, res$n_genes, length(all_genes)))
  all_dfs[[gsm]] <- res$df
}
full_df <- bind_rows(all_dfs)
cat("\n  Combined rows:", nrow(full_df), "\n")

# ---- Merge metadata ----
cat("\n[3/9] Merging with metadata...\n")
full_df <- full_df %>%
  left_join(mol3 %>% select(GSM, barcode_clean, orig.ident, subpop, stage),
            by = c("GSM","barcode_clean")) %>%
  filter(!is.na(subpop))
cat("  Rows after filter:", nrow(full_df), "\n")

# ---- CPM ----
cat("\n[4/9] CPM per gene...\n")
genes_in_df <- intersect(all_genes, colnames(full_df))
cat("  Genes in dataframe:", length(genes_in_df), "/", length(all_genes), "\n")
miss <- setdiff(all_genes, genes_in_df)
if (length(miss) > 0) cat("  Missing genes:", paste(miss, collapse = ", "), "\n")

for (g in genes_in_df) {
  full_df[[paste0(g,"_cpm")]] <- full_df[[g]] / full_df$total_umi * 1e6
  full_df[[paste0(g,"_log")]] <- log1p(full_df[[paste0(g,"_cpm")]])
}

# ---- Module scores ----
cat("\n[5/9] Module scores per cell...\n")
module_names <- character(0)
for (m in names(pathways)) {
  gp <- intersect(pathways[[m]], genes_in_df)
  if (length(gp) < 1) { cat(sprintf("  [SKIP] %-30s: 0 genes\n", m)); next }
  cols <- paste0(gp, "_log")
  full_df[[m]] <- rowMeans(full_df[, cols, drop = FALSE], na.rm = TRUE)
  module_names <- c(module_names, m)
  cat(sprintf("  [OK]   %-30s: %d/%d genes, mean=%.3f\n",
              m, length(gp), length(pathways[[m]]), mean(full_df[[m]], na.rm = TRUE)))
}
stopifnot(length(module_names) > 0)

# ---- Per-mouse aggregation ----
cat("\n[6/9] Per-mouse × subpop...\n")
per_mouse <- full_df %>%
  group_by(orig.ident, subpop) %>%
  summarise(across(all_of(module_names), ~mean(.x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 5)
cat("  Per-mouse rows (n_cells >= 5):", nrow(per_mouse), "\n")

# ---- lmer + paired ----
cat("\n[7/9] lmer + paired t-test per module...\n")
results_list <- list()
for (m in module_names) {
  d_m <- per_mouse[, c("orig.ident","subpop", m)]; names(d_m)[3] <- "score"
  d_m$subpop <- factor(d_m$subpop, levels = c("MOL2_like","MOL56_like"))
  n_mol2  <- sum(d_m$subpop == "MOL2_like")
  n_mol56 <- sum(d_m$subpop == "MOL56_like")

  lmer_res <- tryCatch({
    fit <- lmerTest::lmer(score ~ subpop + (1 | orig.ident), data = d_m, REML = FALSE)
    coefs <- summary(fit)$coefficients
    if ("subpopMOL56_like" %in% rownames(coefs)) {
      list(p = coefs["subpopMOL56_like","Pr(>|t|)"],
           coef = coefs["subpopMOL56_like","Estimate"],
           se = coefs["subpopMOL56_like","Std. Error"],
           note = "lmer_ok")
    } else list(p = NA, coef = NA, se = NA, note = "no_subpop_term")
  }, error = function(e) list(p = NA, coef = NA, se = NA,
                              note = paste0("lmer_failed: ", conditionMessage(e))))

  wide <- d_m %>% pivot_wider(names_from = subpop, values_from = score)
  both <- wide %>% filter(!is.na(MOL2_like) & !is.na(MOL56_like))
  n_paired <- nrow(both)
  paired_res <- if (n_paired >= 3) {
    tryCatch({
      tt <- t.test(both$MOL56_like, both$MOL2_like, paired = TRUE)
      list(p = tt$p.value, coef = mean(both$MOL56_like - both$MOL2_like))
    }, error = function(e) list(p = NA, coef = NA))
  } else list(p = NA, coef = NA)

  results_list[[m]] <- data.frame(
    module = m,
    n_genes_present = length(intersect(pathways[[m]], genes_in_df)),
    n_genes_total   = length(pathways[[m]]),
    n_mice_mol2     = n_mol2,
    n_mice_mol56    = n_mol56,
    n_paired        = n_paired,
    lmer_coef       = lmer_res$coef,
    lmer_se         = lmer_res$se,
    lmer_p          = lmer_res$p,
    lmer_note       = lmer_res$note,
    paired_coef     = paired_res$coef,
    paired_p        = paired_res$p,
    stringsAsFactors = FALSE
  )
}
results <- bind_rows(results_list)
results$lmer_p_adj_BH   <- p.adjust(results$lmer_p,   "BH")
results$paired_p_adj_BH <- p.adjust(results$paired_p, "BH")

cat("\n[8/9] RESULTS:\n")
print(results[, c("module","n_mice_mol2","n_mice_mol56","lmer_coef","lmer_p",
                  "lmer_p_adj_BH","paired_p","paired_p_adj_BH")], row.names = FALSE)

# ---- Comparison with per-cell ----
cat("\n[9/9] Comparison vs per-cell (09_pathway)...\n")
old_cand <- c("results/snrnaseq/09_pathway_forest_data.csv",
              "../results/snrnaseq/09_pathway_forest_data.csv")
old_path <- old_cand[file.exists(old_cand)][1]
if (!is.na(old_path)) {
  old <- fread(old_path, data.table = FALSE)[, c("module","coef_MOL56_vs_MOL2","pval")]
  names(old) <- c("module","percell_coef","percell_p")
  comparison <- merge(results[, c("module","lmer_coef","lmer_p","paired_p")],
                      old, by = "module", all.x = TRUE)
  comparison$delta_coef <- comparison$lmer_coef - comparison$percell_coef
  print(comparison, row.names = FALSE)
  write.csv(comparison, "results/snrnaseq/09b_percell_vs_permouse.csv", row.names = FALSE)
  cat("[OK] results/snrnaseq/09b_percell_vs_permouse.csv\n")
}

write.csv(results,   "results/snrnaseq/09b_pathway_per_mouse_pvals.csv", row.names = FALSE)
write.csv(per_mouse, "results/snrnaseq/09b_per_mouse_module_scores.csv", row.names = FALSE)
cat("\n[OK] results/snrnaseq/09b_pathway_per_mouse_pvals.csv\n")
cat("[OK] results/snrnaseq/09b_per_mouse_module_scores.csv\n")
cat("\n=== FIN ===\n")
