# ============================================================
# 32_verify_liver_numbers.R
# Verifica log2FC de SMPD3 y SMS1 contra WT/KO reportados
# Uso: Rscript scripts/32_verify_liver_numbers.R
# Salida: logs/32_verify_liver_log.txt, results/32_liver_check.csv
# ============================================================

suppressPackageStartupMessages({ library(dplyr) })
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/32_verify_liver_log.txt", split = TRUE)
on.exit(sink())

cat("=== VERIFICACION DE GENES HEPATICOS ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# --- Localizar CSV DESeq2 de SRP074136 ---
f <- c("results/SRP074136_DESeq2_inflamed_vs_WT.csv",
       "results/SRP074136_DESeq2_inflamed_vs_WT_annotated.csv")
f <- f[file.exists(f)][1]
if (is.na(f)) {
  stop("No encuentro el CSV DESeq2 de SRP074136. Verifica results/")
}
cat("Fuente:", f, "\n\n")

deseq <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
cat("Columnas:", paste(colnames(deseq), collapse = ", "), "\n\n")

# Reemplaza la sección "detectar columna de gen" en 32_verify_liver_numbers.R
# Añadir mapping Ensembl → Symbol:
suppressPackageStartupMessages(library(org.Mm.eg.db))
deseq$symbol <- mapIds(org.Mm.eg.db,
                       keys = deseq$gene_id,
                       keytype = "ENSEMBL",
                       column = "SYMBOL")
gene_col <- "symbol"

# --- Genes de interes ---
genes <- c("Smpd3", "Sgms1", "Sgms2", "B4galt5", "B4galt6", "Ugcg", "Cers2")
sub <- deseq[deseq[[gene_col]] %in% genes, ]
print(sub[, intersect(c(gene_col, "log2FoldChange", "baseMean", "padj"),
                      colnames(sub))])
cat("\n")

# --- Recalcular log2FC de WT/KO reportados en la tabla del paper ---
paper_tbl <- data.frame(
  gene     = c("Smpd3", "Sgms1", "Sgms2", "B4galt5", "B4galt6", "Cers2"),
  WT       = c(0.8,   17,    61,    27,    1.1,   159),
  KO       = c(22.1,  25,    70,    93,    6.0,   159),
  log2FC_p = c(3.68,  0.51,  0.21,  1.74,  1.72,  -0.003)
)
paper_tbl$log2FC_calc <- log2(paper_tbl$KO / paper_tbl$WT)
paper_tbl$delta      <- paper_tbl$log2FC_calc - paper_tbl$log2FC_p

cat("=== VERIFICACION TABLA DEL PAPER (WT/KO) ===\n")
print(paper_tbl)
cat("\n")

# --- Cruce con DESeq2 ---
cat("=== COMPARACION CON DESeq2 ===\n")
for (i in seq_len(nrow(paper_tbl))) {
  g <- paper_tbl$gene[i]
  row <- sub[sub[[gene_col]] == g, ]
  if (nrow(row) == 0) {
    cat(sprintf("  %-8s: NO en DESeq2\n", g))
    next
  }
  lfc_deseq <- row$log2FoldChange[1]
  lfc_paper <- paper_tbl$log2FC_p[i]
  lfc_wtko  <- paper_tbl$log2FC_calc[i]
  cat(sprintf("  %-8s  paper=%+.3f  WT/KO=%+.3f  DESeq2=%+.3f  |  paper vs WT/KO=%+.3f  paper vs DESeq2=%+.3f\n",
              g, lfc_paper, lfc_wtko, lfc_deseq,
              lfc_paper - lfc_wtko, lfc_paper - lfc_deseq))
}

cat("\n=== INTERPRETACION ===\n")
cat("  Si 'paper vs WT/KO' > 1.0 -> la tabla del paper tiene valores WT/KO\n")
cat("  que no coinciden con el log2FC reportado.\n")
cat("  Si 'paper vs DESeq2' > 0.5 -> el log2FC del paper NO coincide con DESeq2.\n")
cat("  En ese caso, usar DESeq2 (que es la fuente primaria) en la tabla final.\n")

# --- Guardar ---
out <- merge(paper_tbl, sub[, c(gene_col, "log2FoldChange", "baseMean", "padj")],
             by.x = "gene", by.y = gene_col, all.x = TRUE)
write.csv(out, "results/32_liver_check.csv", row.names = FALSE)
cat("\n  -> results/32_liver_check.csv\n")
cat("\n=== FIN ===\n")
