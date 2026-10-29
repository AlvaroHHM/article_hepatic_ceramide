# ============================================================
# 25_reanalysis_GSE253318_microglia.R
# Reanalisis de GSE253318: microglia EAE vs control (FPKM)
# Objetivo: comparar la via de esfingolipidos en microglia
# con astrocitos (GSE139531) y bulk SNC (SRP071304)
# Uso: Rscript scripts/25_reanalysis_GSE253318_microglia.R
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(tidyr)
})

dir.create("data/GSE253318", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/25_GSE253318_log.txt", split = TRUE)
on.exit(sink())

message("=== GSE253318: microglia EAE vs control ===")
message("Fecha: ", Sys.time())

# ------------------------------------------------------------
# 1. Localizar archivo FPKM
# ------------------------------------------------------------
supp_dir <- "data/GSE253318"
if (!dir.exists(supp_dir) || length(list.files(supp_dir, recursive = TRUE)) == 0) {
  message("Descargando archivos suplementarios de GSE253318...")
  getGEOSuppFiles("GSE253318", baseDir = supp_dir, makeDirectory = TRUE)
}
supp_files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE)

# Preferir FPKM sobre counts
fpkm_candidates <- supp_files[grepl("FPKM", basename(supp_files), ignore.case = TRUE)]
if (length(fpkm_candidates) == 0) {
  fpkm_candidates <- supp_files[grepl("\\.(txt|tsv|csv)(\\.gz)?$", basename(supp_files), ignore.case = TRUE)]
}
fpkm_file <- fpkm_candidates[1]
message("Archivo: ", basename(fpkm_file))

# ------------------------------------------------------------
# 2. Leer la matriz de expresion
# ------------------------------------------------------------
if (grepl("\\.gz$", fpkm_file)) {
  con <- gzfile(fpkm_file, "r")
  expr_raw <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
  close(con)
} else {
  expr_raw <- read.delim(fpkm_file, stringsAsFactors = FALSE, check.names = FALSE)
}
message("Dimensiones: ", paste(dim(expr_raw), collapse = " x "))

# Primera columna = gen
rownames(expr_raw) <- expr_raw[[1]]
expr_raw <- expr_raw[, -1, drop = FALSE]
expr_raw[] <- lapply(expr_raw, function(x) as.numeric(as.character(x)))

# ------------------------------------------------------------
# 3. Definir grupos (case = EAE, con = Control)
# ------------------------------------------------------------
group <- ifelse(grepl("^case", colnames(expr_raw), ignore.case = TRUE), "EAE",
         ifelse(grepl("^con",  colnames(expr_raw), ignore.case = TRUE), "Control",
                NA_character_))
message("\nGrupos detectados:")
print(table(group, useNA = "ifany"))

valid <- !is.na(group)
expr_raw <- expr_raw[, valid, drop = FALSE]
group <- group[valid]

if (sum(group == "EAE") < 2 || sum(group == "Control") < 2) {
  stop("Menos de 2 muestras por grupo. Revisar nombres de columna.")
}
message("\nMuestras EAE: ", sum(group == "EAE"),
        " | Muestras Control: ", sum(group == "Control"))

# ------------------------------------------------------------
# 4. Genes de interes
# ------------------------------------------------------------
genes_interes <- c(
  "B4galt6", "B4galt5", "Ugcg", "Smpd1", "Smpd3", "Smpd2", "Smpd4",
  "Sgms1", "Sgms2", "Cers1", "Cers2", "Cers3", "Cers4", "Cers5", "Cers6",
  "Sptlc1", "Sptlc2", "Kdsr", "Degs1", "Degs2",
  "Asah1", "Asah2", "Acer1", "Acer2", "Acer3",
  "Sphk1", "Sphk2", "Sgpl1", "S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5",
  "Ormdl1", "Ormdl2", "Ormdl3", "Fa2h",
  "Pla2g4a", "Mavs", "Nfkb1", "Rela", "Irf3", "Irf7"
)

genes_presentes <- intersect(genes_interes, rownames(expr_raw))
message("\nGenes de interes presentes: ", length(genes_presentes), "/", length(genes_interes))
message("Genes presentes: ", paste(genes_presentes, collapse = ", "))
message("Genes ausentes: ", paste(setdiff(genes_interes, genes_presentes), collapse = ", "))

expr_sub <- expr_raw[genes_presentes, , drop = FALSE]

# ------------------------------------------------------------
# 5. Medias por grupo y log2FC
# ------------------------------------------------------------
res <- data.frame(gene = genes_presentes, stringsAsFactors = FALSE)
for (g in c("Control", "EAE")) {
  res[[paste0("mean_", g)]] <- apply(expr_sub, 1, function(x) mean(x[group == g], na.rm = TRUE))
}
res$log2FC <- log2(res$mean_EAE + 1) - log2(res$mean_Control + 1)
res <- res[order(-res$log2FC), ]

write.csv(res, "results/25_GSE253318_microglia_genes.csv", row.names = FALSE)
message("OK results/25_GSE253318_microglia_genes.csv")

message("\n========================================")
message("RPKM MEDIO POR GEN Y GRUPO (microglia)")
message("========================================")
print(res, row.names = FALSE, digits = 2)

# ------------------------------------------------------------
# 6. Interpretacion para el nodo 7
# ------------------------------------------------------------
message("\n========================================")
message("INTERPRETACION NODO 7 (MICROGLIA)")
message("========================================")

b4 <- res$log2FC[res$gene == "B4galt6"]
if (length(b4) > 0) {
  message(sprintf("\nB4galt6 log2FC (EAE vs Control): %+.2f", b4))
  if (b4 > 1) message("  -> B4galt6 se induce en microglia EAE")
  else if (b4 < -1) message("  -> B4galt6 se suprime en microglia EAE")
  else message("  -> B4galt6 no cambia significativamente en microglia EAE")
}

smpd3 <- res$log2FC[res$gene == "Smpd3"]
sgms1 <- res$log2FC[res$gene == "Sgms1"]
if (length(smpd3) > 0 && length(sgms1) > 0) {
  message(sprintf("\nSmpd3 log2FC: %+.2f | Sgms1 log2FC: %+.2f", smpd3, sgms1))
  if (smpd3 > 1 && sgms1 < 1) {
    message("  -> Imbalance Smpd3 > Sgms1 en microglia (similar al higado)")
  } else {
    message("  -> No hay imbalance Smpd3 > Sgms1 en microglia")
  }
}

message("\nVia S1P en microglia:")
for (g in c("Asah1", "Sgpl1", "S1pr1", "S1pr3", "Sphk1", "Sphk2")) {
  fc <- res$log2FC[res$gene == g]
  if (length(fc) > 0) message(sprintf("  %s: log2FC = %+.2f", g, fc))
}

message("\n=== FIN ===")
message("Log: logs/25_GSE253318_log.txt")
