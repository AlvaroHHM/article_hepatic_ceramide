# ============================================================
# 04_verify_liver_numbers.R
# Extracts log2FC and padj for sphingolipid cycle genes in
# SRP074136 (Mdr2-KO liver) from the DESeq2 output.
#
# Input:
#   results/bulk/SRP074136_DESeq2_inflamed_vs_WT.csv
# Output:
#   results/bulk/32_liver_check.csv
#   logs/32_verify_liver_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

dir.create("results", showWarnings = FALSE)
dir.create("results/bulk", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE)

sink("logs/32_verify_liver_log.txt", split = TRUE)
on.exit(sink())

cat("=== LIVER GENE VALUES (SRP074136, DESeq2) ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# --- Locate DESeq2 output ---
f <- c("results/bulk/SRP074136_DESeq2_inflamed_vs_WT.csv",
       "results/SRP074136_DESeq2_inflamed_vs_WT.csv")
f <- f[file.exists(f)][1]
if (is.na(f)) {
  stop("DESeq2 output not found. Expected results/bulk/SRP074136_DESeq2_inflamed_vs_WT.csv")
}
cat("Source:", f, "\n\n")

deseq <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
cat("Columns:", paste(colnames(deseq), collapse = ", "), "\n\n")

# --- Map Ensembl (no version) -> Symbol ---
deseq$symbol <- suppressMessages(
  AnnotationDbi::mapIds(org.Mm.eg.db,
                        keys = deseq$ensembl_clean,
                        keytype = "ENSEMBL",
                        column = "SYMBOL",
                        multiVals = "first")
)
cat("Mapped:", sum(!is.na(deseq$symbol)), "/", nrow(deseq), "\n\n")

# --- Genes of interest (explicit dplyr:: calls to avoid AnnotationDbi::select mask) ---
genes <- c("Smpd3", "Sgms1", "Sgms2", "B4galt5", "B4galt6", "Ugcg", "Cers2")

sub <- deseq %>%
  dplyr::filter(symbol %in% genes) %>%
  dplyr::select(gene = symbol, log2FoldChange, baseMean, padj) %>%
  dplyr::arrange(gene)

cat("=== Sphingolipid cycle genes in SRP074136 (DESeq2) ===\n")
print(sub, digits = 4)
cat("\n")

# --- Save ---
write.csv(sub, "results/bulk/32_liver_check.csv", row.names = FALSE)
cat("-> results/bulk/32_liver_check.csv\n")
cat("\n=== FIN ===\n")
