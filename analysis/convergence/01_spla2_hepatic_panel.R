# ============================================================
# 35_spla2_hepatic_panel.R
# Panel extendido de sPLA2 secretoras en higado EAE (GSE312673)
# + reanalisis del outlier de colon (Takahama)
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

sink("logs/35_spla2_hepatic_log.txt", split = TRUE)
on.exit(sink())

cat("============================================================\n")
cat("sPLA2 SECRETED PANEL IN LIVER EAE + COLON OUTLIER CHECK\n")
cat("Date:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# ============================================================
# PARTE A: PANEL DE sPLA2 EN HIGADO EAE (GSE312673)
# ============================================================
cat("=== PARTE A: sPLA2 secretadas en higado EAE ===\n\n")

count_file <- "data/GSE312673/GSE312673_counts_matrix.tsv.gz"
if (!file.exists(count_file)) stop("Missing: ", count_file)

counts_raw <- read.delim(count_file, row.names = 1, check.names = FALSE,
                          stringsAsFactors = FALSE, comment.char = "#")
cat("Counts:", nrow(counts_raw), "x", ncol(counts_raw), "\n")

ids_clean <- sub("\\..*$", "", rownames(counts_raw))
symbols <- mapIds(org.Mm.eg.db,
                  keys = ids_clean, keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
counts_raw$symbol <- as.character(symbols)
counts_raw <- counts_raw[!is.na(counts_raw$symbol) &
                          nzchar(counts_raw$symbol), ]

samples <- setdiff(colnames(counts_raw), "symbol")
counts_sym <- counts_raw %>%
  group_by(symbol) %>%
  summarise(across(all_of(samples), ~sum(.x, na.rm = TRUE)), .groups = "drop")
counts_mat <- as.data.frame(counts_sym)
rownames(counts_mat) <- counts_mat$symbol
counts_mat$symbol <- NULL
cat("Genes agregados:", nrow(counts_mat), "\n")

parse_group <- function(s) {
  case_when(
    grepl("peak.*EAE|EAE.*peak", s, ignore.case = TRUE) ~ "Peak_EAE",
    grepl("peak.*NC|NC.*peak",   s, ignore.case = TRUE) ~ "Peak_Ctrl",
    grepl("pre.*EAE|EAE.*pre",   s, ignore.case = TRUE) ~ "Pre_EAE",
    grepl("pre.*NC|NC.*pre",     s, ignore.case = TRUE) ~ "Pre_Ctrl",
    TRUE ~ NA_character_
  )
}
colData <- data.frame(
  sample = colnames(counts_mat),
  group  = factor(parse_group(colnames(counts_mat)),
                  levels = c("Pre_Ctrl", "Pre_EAE", "Peak_Ctrl", "Peak_EAE"))
)
rownames(colData) <- colData$sample
colData$sample <- NULL
print(table(colData$group))

keep <- rowSums(counts_mat) >= 10
counts_filt <- counts_mat[keep, ]
cat("Genes filtrados:", sum(keep), "/", nrow(counts_mat), "\n")

dds <- DESeqDataSetFromMatrix(counts_filt, colData, design = ~ group)
dds <- DESeq(dds, quiet = TRUE)

# --- Panel extendido de sPLA2 secretoras ---
spla2_panel <- c(
  # Grupo IB (pancreatica, digestiva)
  "Pla2g1b",
  # Grupo IIA (reactante de fase aguda, AMP principal)
  "Pla2g2a",
  # Grupo IID/IIE/IIF (inflamatorias, mucosa)
  "Pla2g2d", "Pla2g2e", "Pla2g2f",
  # Grupo V (inflamatoria, sepsis)
  "Pla2g5",
  # Grupo III, X, XIIA
  "Pla2g3", "Pla2g10", "Pla2g12a", "Pla2g12b"
)

# Extraer resultados
contrasts <- list(
  Peak_EAE_vs_Ctrl = c("group", "Peak_EAE", "Peak_Ctrl"),
  Pre_EAE_vs_Ctrl  = c("group", "Pre_EAE",  "Pre_Ctrl")
)

spla2_results <- list()
for (nm in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[nm]], alpha = 0.05)
  df  <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[, c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")]
  df$sig <- ifelse(is.na(df$padj), "NA",
                   ifelse(df$padj < 0.001, "***",
                          ifelse(df$padj < 0.01, "**",
                                 ifelse(df$padj < 0.05, "*", "ns"))))
  spla2_results[[nm]] <- df %>%
    dplyr::filter(gene %in% spla2_panel) %>%
    dplyr::arrange(desc(abs(log2FoldChange)))
}

cat("\n=== sPLA2 PANEL: Peak EAE vs Ctrl ===\n")
print(as.data.frame(spla2_results$Peak_EAE_vs_Ctrl))
cat("\n=== sPLA2 PANEL: Pre EAE vs Ctrl ===\n")
print(as.data.frame(spla2_results$Pre_EAE_vs_Ctrl))

# Guardar
write.csv(spla2_results$Peak_EAE_vs_Ctrl,
          "results/convergence/35_sPLA2_peak_EAE_vs_ctrl.csv", row.names = FALSE)
write.csv(spla2_results$Pre_EAE_vs_Ctrl,
          "results/convergence/35_sPLA2_pre_EAE_vs_ctrl.csv", row.names = FALSE)

# Veredicto sobre Pla2g2a (el AMP real)
cat("\n=== VEREDICTO: Pla2g2a (sPLA2-IIA, el AMP principal) ===\n")
p2g2a <- spla2_results$Peak_EAE_vs_Ctrl %>%
  dplyr::filter(gene == "Pla2g2a")
if (nrow(p2g2a) > 0) {
  cat(sprintf("  Peak: log2FC = %+.2f, padj = %.3g (%s)\n",
              p2g2a$log2FoldChange[1], p2g2a$padj[1], p2g2a$sig[1]))
  if (!is.na(p2g2a$padj[1]) && p2g2a$padj[1] < 0.05 &&
      p2g2a$log2FoldChange[1] > 1) {
    cat("  -> Pla2g2a SE INDUCE en higado EAE (respuesta inflamatoria)\n")
    cat("     Esto es consistente con fase aguda. El AMP SI se activa\n")
    cat("     en higado, aunque la digestiva (Pla2g1b) no se compense.\n")
  } else if (!is.na(p2g2a$padj[1]) && p2g2a$padj[1] < 0.05 &&
             p2g2a$log2FoldChange[1] < -1) {
    cat("  -> Pla2g2a SE SUPRIME en higado EAE. El AMP hepatico cae.\n")
  } else {
    cat("  -> Pla2g2a NO cambia significativamente en higado EAE.\n")
  }
}

# Veredicto sobre Pla2g1b
cat("\n=== VEREDICTO: Pla2g1b (sPLA2-IB, digestiva) ===\n")
p2g1b <- spla2_results$Peak_EAE_vs_Ctrl %>%
  dplyr::filter(gene == "Pla2g1b")
if (nrow(p2g1b) > 0) {
  cat(sprintf("  Peak: log2FC = %+.2f, padj = %.3g (%s)\n",
              p2g1b$log2FoldChange[1], p2g1b$padj[1], p2g1b$sig[1]))
  if (!is.na(p2g1b$log2FoldChange[1])) {
    if (p2g1b$log2FoldChange[1] > 1 &&
        !is.na(p2g1b$padj[1]) && p2g1b$padj[1] < 0.05) {
      cat("  -> Pla2g1b SE INDUCE en higado EAE (expresion ectopica)\n")
      cat("     Seria una respuesta compensatoria del higado por la\n")
      cat("     perdida pancreatica/intestinal. NO esperado a priori.\n")
    } else {
      cat("  -> Pla2g1b NO se induce en higado EAE (esperado;\n")
      cat("     es enzima pancreatica/intestinal, no hepatica)\n")
    }
  }
}

cat("\n=== OUTPUTS ===\n")
cat("  results/convergence/35_sPLA2_peak_EAE_vs_ctrl.csv\n")
cat("  results/convergence/35_sPLA2_pre_EAE_vs_ctrl.csv\n")
cat("  logs/35_spla2_hepatic_log.txt\n")
cat("\n=== FIN ===\n")
