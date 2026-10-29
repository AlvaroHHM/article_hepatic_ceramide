# ============================================================
# 34_DESeq2_liver_EAE_GSE312673.R
# DESeq2 analysis of liver transcriptome in EAE (GSE312673)
# Contrasts: Peak_EAE vs Peak_Ctrl, Pre_EAE vs Pre_Ctrl
# Requires: data/GSE312673/GSE312673_counts_matrix.tsv.gz
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

sink("logs/34_DESeq2_liver_EAE_log.txt", split = TRUE)
on.exit(sink())

cat("============================================================\n")
cat("DESEQ2 ANALYSIS OF LIVER IN EAE (GSE312673)\n")
cat("Date:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load counts
# ------------------------------------------------------------
count_file <- "data/GSE312673/GSE312673_counts_matrix.tsv.gz"
if (!file.exists(count_file)) stop("Missing: ", count_file)

counts_raw <- read.delim(count_file, row.names = 1, check.names = FALSE,
                          stringsAsFactors = FALSE, comment.char = "#")
cat("Counts matrix:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("Samples:", paste(colnames(counts_raw), collapse = ", "), "\n\n")

# ------------------------------------------------------------
# 2. Map Ensembl -> Symbol and aggregate
# ------------------------------------------------------------
ids_clean <- sub("\\..*$", "", rownames(counts_raw))
symbols <- mapIds(org.Mm.eg.db,
                  keys      = ids_clean,
                  keytype   = "ENSEMBL",
                  column    = "SYMBOL",
                  multiVals = "first")
cat("Symbols mapped:", sum(!is.na(symbols)), "/", length(symbols), "\n")

counts_raw$symbol <- as.character(symbols)
counts_raw <- counts_raw[!is.na(counts_raw$symbol) & nzchar(counts_raw$symbol), ]

samples <- setdiff(colnames(counts_raw), "symbol")
counts_sym <- counts_raw %>%
  group_by(symbol) %>%
  summarise(across(all_of(samples), ~sum(.x, na.rm = TRUE)), .groups = "drop")
counts_mat <- as.data.frame(counts_sym)
rownames(counts_mat) <- counts_mat$symbol
counts_mat$symbol <- NULL
cat("Aggregated genes:", nrow(counts_mat), "\n\n")

# ------------------------------------------------------------
# 3. Build colData
# ------------------------------------------------------------
parse_group <- function(s) {
  case_when(
    grepl("peak.*EAE|EAE.*peak", s, ignore.case = TRUE) ~ "Peak_EAE",
    grepl("peak.*NC|NC.*peak",   s, ignore.case = TRUE) ~ "Peak_Ctrl",
    grepl("pre.*EAE|EAE.*pre",   s, ignore.case = TRUE) ~ "Pre_EAE",
    grepl("pre.*NC|NC.*pre",     s, ignore.case = TRUE) ~ "Pre_Ctrl",
    TRUE ~ NA_character_
  )
}
sample_groups <- parse_group(colnames(counts_mat))
cat("Group assignment:\n")
print(table(sample_groups))
if (any(is.na(sample_groups))) {
  stop("Unassigned samples: ",
       paste(colnames(counts_mat)[is.na(sample_groups)], collapse = ", "))
}

colData <- data.frame(
  sample = colnames(counts_mat),
  group  = factor(sample_groups,
                  levels = c("Pre_Ctrl", "Pre_EAE",
                             "Peak_Ctrl", "Peak_EAE"))
)
rownames(colData) <- colData$sample
colData$sample <- NULL

# ------------------------------------------------------------
# 4. Filter low-count genes
# ------------------------------------------------------------
keep <- rowSums(counts_mat) >= 10
cat("\nGenes passing rowSums >= 10:", sum(keep), "/", nrow(counts_mat), "\n")
counts_filt <- counts_mat[keep, ]

# ------------------------------------------------------------
# 5. Run DESeq2
# ------------------------------------------------------------
cat("\n--- Running DESeq2 ---\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = colData,
  design    = ~ group
)
dds <- DESeq(dds, quiet = TRUE)
cat("resultsNames(dds):\n")
print(resultsNames(dds))
cat("\n")

# ------------------------------------------------------------
# 6. Extract contrasts
# ------------------------------------------------------------
contrasts <- list(
  Peak_EAE_vs_Ctrl = c("group", "Peak_EAE", "Peak_Ctrl"),
  Pre_EAE_vs_Ctrl  = c("group", "Pre_EAE",  "Pre_Ctrl")
)

results_list <- list()
for (nm in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[nm]], alpha = 0.05)
  df  <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[, c("gene", "baseMean", "log2FoldChange",
               "lfcSE", "stat", "pvalue", "padj")]
  df$sig <- ifelse(is.na(df$padj), "NA",
                   ifelse(df$padj < 0.001, "***",
                          ifelse(df$padj < 0.01, "**",
                                 ifelse(df$padj < 0.05, "*", "ns"))))
  results_list[[nm]] <- df
  out_file <- paste0("results/34_DESeq2_liver_EAE_", nm, ".csv")
  write.csv(df, out_file, row.names = FALSE)
  n_sig <- sum(df$sig %in% c("*", "**", "***"), na.rm = TRUE)
  cat(sprintf("  %-20s: %d genes, %d significant (padj < 0.05) -> %s\n",
              nm, nrow(df), n_sig, out_file))
}

# ------------------------------------------------------------
# 7. Panel of interest
# ------------------------------------------------------------
panel_genes <- c(
  "Smpd3", "Sgms1", "Sgms2",
  "B4galt5", "B4galt6", "Ugcg",
  "Cers2", "Cers5", "Cers6",
  "Cyp7a1", "Cyp8b1", "Nr1h4", "Abcb11", "Slc10a1",
  "Tnf", "Il6", "Il1b", "Ccl2",
  "Sphk1", "S1pr1"
)

panel_df <- tibble(gene = panel_genes)
for (nm in names(results_list)) {
  df  <- results_list[[nm]]
  sel <- df[df$gene %in% panel_genes, c("gene", "log2FoldChange",
                                         "padj", "sig", "baseMean")]
  colnames(sel) <- c("gene",
                     paste0("log2FC_", nm),
                     paste0("padj_",    nm),
                     paste0("sig_",     nm),
                     paste0("baseMean_", nm))
  panel_df <- left_join(panel_df, sel, by = "gene")
}

# Which panel genes were filtered out (expression < 10)?
panel_missing <- setdiff(panel_genes, rownames(counts_filt))
if (length(panel_missing) > 0) {
  cat("\nPanel genes filtered out (rowSums < 10):",
      paste(panel_missing, collapse = ", "), "\n")
}

cat("\n=== PANEL OF INTEREST ===\n")
print(as.data.frame(panel_df))
write.csv(panel_df, "results/34_DESeq2_liver_EAE_panel.csv", row.names = FALSE)

# ------------------------------------------------------------
# 8. Side-by-side with Mdr2-KO
# ------------------------------------------------------------
mdr2_reference <- tribble(
  ~gene,    ~log2FC_Mdr2KO, ~padj_Mdr2KO,
  "Smpd3",  4.80,           2e-58,
  "Sgms1",  0.01,           0.98,
  "Sgms2",  -0.23,          0.15,
  "B4galt6",3.77,           1e-43,
  "B4galt5",1.89,           1e-20,
  "Cers2",  -0.15,          0.43
)

comparison <- mdr2_reference %>%
  left_join(panel_df %>%
              select(gene,
                     log2FC_Peak_EAE_vs_Ctrl,
                     padj_Peak_EAE_vs_Ctrl,
                     log2FC_Pre_EAE_vs_Ctrl,
                     padj_Pre_EAE_vs_Ctrl),
            by = "gene")

cat("\n=== MDR2-KO (colestasis) VS EAE (enfermedad) ===\n")
cat(sprintf("%-10s %12s %10s %12s %10s %12s %10s\n",
            "gene", "Mdr2KO_lfc", "Mdr2KO_p", "EAE_Peak_lfc",
            "EAE_Peak_p", "EAE_Pre_lfc", "EAE_Pre_p"))
for (i in seq_len(nrow(comparison))) {
  row <- comparison[i, ]
  cat(sprintf("%-10s %+12.2f %10.2g %+12.2f %10.2g %+12.2f %10.2g\n",
              row$gene,
              row$log2FC_Mdr2KO, row$padj_Mdr2KO,
              row$log2FC_Peak_EAE_vs_Ctrl, row$padj_Peak_EAE_vs_Ctrl,
              row$log2FC_Pre_EAE_vs_Ctrl,  row$padj_Pre_EAE_vs_Ctrl))
}
write.csv(comparison, "results/34_comparison_Mdr2KO_vs_EAE.csv",
          row.names = FALSE)

# ------------------------------------------------------------
# 9. Interpretation
# ------------------------------------------------------------
cat("\n=== INTERPRETATION ===\n")
check_gene <- function(g) {
  if (!(g %in% panel_df$gene)) return(NULL)
  row <- panel_df[panel_df$gene == g, ]
  list(
    pre_lfc  = row$log2FC_Pre_EAE_vs_Ctrl,
    pre_padj = row$padj_Pre_EAE_vs_Ctrl,
    pre_sig  = row$sig_Pre_EAE_vs_Ctrl,
    pk_lfc   = row$log2FC_Peak_EAE_vs_Ctrl,
    pk_padj  = row$padj_Peak_EAE_vs_Ctrl,
    pk_sig   = row$sig_Peak_EAE_vs_Ctrl
  )
}

for (g in c("Smpd3", "Sgms1", "B4galt6", "Cyp7a1", "Cyp8b1", "Il1b", "Ccl2")) {
  r <- check_gene(g)
  if (is.null(r)) next
  cat(sprintf("  %-8s  Pre: lfc=%+5.2f (%s)  Peak: lfc=%+5.2f (%s)\n",
              g,
              r$pre_lfc, r$pre_sig,
              r$pk_lfc,  r$pk_sig))
}

# Final verdict
cat("\n=== VERDICT ===\n")
smpd3 <- check_gene("Smpd3")
sgms1 <- check_gene("Sgms1")
if (!is.null(smpd3) && !is.null(sgms1)) {
  up_smpd3 <- (smpd3$pk_lfc > 0.5 && !is.na(smpd3$pk_padj) && smpd3$pk_padj < 0.1)
  flat_sgms1 <- (abs(sgms1$pk_lfc) < 0.5 ||
                 (is.na(sgms1$pk_padj) || sgms1$pk_padj > 0.1))
  if (up_smpd3 && flat_sgms1) {
    cat("  PATRON CONFIRMADO en EAE: SMPD3 inducido, SMS1 plano.\n")
    cat("  El imbalance hepático productor/procesador se replica en el\n")
    cat("  modelo de enfermedad (EAE), no solo en colestasis (Mdr2-KO).\n")
  } else if (up_smpd3 && !flat_sgms1) {
    cat("  PARCIAL: SMPD3 inducido, pero SMS1 también se mueve.\n")
    cat("  El imbalance es menos selectivo en EAE que en Mdr2-KO.\n")
  } else {
    cat("  NO CONFIRMADO: el patrón hepático difiere entre EAE y Mdr2-KO.\n")
  }
}

# Check biliar
cyp <- check_gene("Cyp7a1")
if (!is.null(cyp)) {
  cat(sprintf("\n  Cyp7a1 (síntesis biliar) en Peak: lfc=%+.2f (%s)\n",
              cyp$pk_lfc, cyp$pk_sig))
  if (cyp$pk_lfc < -1 && !is.na(cyp$pk_padj) && cyp$pk_padj < 0.1) {
    cat("  -> Supresión biliar CONFIRMADA en hígado EAE.\n")
    cat("     Conecta nodo 4 (supresión biliar) con nodo 5 (ciclo SM).\n")
  }
}

cat("\n=== OUTPUTS ===\n")
cat("  results/34_DESeq2_liver_EAE_Peak_EAE_vs_Ctrl.csv\n")
cat("  results/34_DESeq2_liver_EAE_Pre_EAE_vs_Ctrl.csv\n")
cat("  results/34_DESeq2_liver_EAE_panel.csv\n")
cat("  results/34_comparison_Mdr2KO_vs_EAE.csv\n")
cat("  logs/34_DESeq2_liver_EAE_log.txt\n")
cat("\n=== FIN ===\n")
