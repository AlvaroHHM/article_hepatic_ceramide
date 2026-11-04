# ============================================================
# 01_oligodendrocytes_GSE250589_annotation.R
# Descarga y verificacion de la anotacion de GSE250589
# (Zheng et al. 2025, snRNA-seq multiome de EAE en medula espinal).
#
# NOTA: la matriz de expresion (44.7 GB) NO se descarga (no cabe
# en el hardware actual). Los scripts 05-10 (analisis MOL3_EAE)
# requieren los 28 H5 files, que se descargan con 00_download_h5.sh.
#
# Uso: Rscript analysis/snrnaseq/01_oligodendrocytes_GSE250589_annotation.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

dir.create("logs/snrnaseq", showWarnings = FALSE, recursive = TRUE)
sink("logs/snrnaseq/01_GSE250589_annotation_log.txt", split = TRUE)
on.exit(sink())

SEP <- "============================================================"
cat(SEP, "\n")
cat("GSE250589 - Descarga y verificacion de la anotacion\n")
cat("Fecha:", format(Sys.time()), "\n")
cat(SEP, "\n\n")

# --- 1. Localizar o descargar la anotacion ---
candidates <- c(
  "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz",
  "../data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
)
ann_file <- candidates[file.exists(candidates)][1]

if (is.na(ann_file)) {
  cat("Anotacion no encontrada localmente. Descargando (139 MB)...\n")
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250589/suppl/GSE250589_EAE_multiomics_annotation.csv.gz"
  ann_file <- "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
  dir.create(dirname(ann_file), showWarnings = FALSE, recursive = TRUE)
  download.file(url, ann_file, mode = "wb", quiet = FALSE)
}
cat("Archivo:", ann_file, "(",
    round(file.size(ann_file)/1024/1024, 1), "MB)\n\n")

# --- 2. Cargar ---
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Instala data.table: install.packages('data.table')")
}
suppressPackageStartupMessages(library(data.table))

cat("Leyendo CSV comprimido...\n")
ann <- fread(ann_file, data.table = FALSE)
cat("Dimensiones:", nrow(ann), "x", ncol(ann), "\n\n")

# --- 3. Verificacion minima ---
cat("Columnas clave (verificacion):\n")
for (col in c("cellType", "subcellType", "cellType_merge",
              "model", "orig.ident", "batch")) {
  if (col %in% colnames(ann)) {
    n_unique <- length(unique(ann[[col]]))
    cat(sprintf("  %-15s presente (%d valores unicos)\n", col, n_unique))
  } else {
    cat(sprintf("  %-15s AUSENTE\n", col))
  }
}

cat("\nDistribucion de 'model':\n")
print(table(ann$model, useNA = "ifany"))

cat("\nDistribucion de 'orig.ident' (top 10):\n")
print(head(sort(table(ann$orig.ident), decreasing = TRUE), 10))

# --- 4. Aviso para scripts downstream ---
n_h5 <- length(list.files("data/GSE250589/h5", pattern = "\\.h5$"))
cat("\n", SEP, "\n")
cat("Estado de los H5 files (necesarios para scripts 05-10):\n")
cat("  Encontrados:", n_h5, "/ 28\n")
if (n_h5 < 28) {
  cat("  Ejecutar: bash analysis/snrnaseq/00_download_h5.sh\n")
  cat("  (descarga ~2.3 GB desde GEO, 10-30 min)\n")
} else {
  cat("  Listos para correr scripts 05-10\n")
}

cat("\n", SEP, "\nFIN\n", SEP, "\n")
