# ============================================================
# 29_oligodendrocytes_GSE250589_annotation.R
# Fase 1 SOLO: descarga y verifica la anotacion de GSE250589
# El RDS con la matriz de expresion (44.7 GB) NO se descarga:
# no cabe en el hardware actual (30 GB libres, 15 GB RAM).
# El objetivo es caracterizar la estructura del dataset para
# el anteproyecto y citarlo como "future work" con datos propios.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

dir.create("data/GSE250589", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/29_GSE250589_annotation_log.txt", split = TRUE)
on.exit(sink())

SEP <- "============================================================"
cat(SEP, "\n")
cat("GSE250589 - Fase 1: verificacion de la anotacion\n")
cat("Fecha:", format(Sys.time()), "\n")
cat(SEP, "\n\n")

# --- 1. Descargar solo si falta ---
ann_file <- "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
if (!file.exists(ann_file)) {
  cat("Descargando anotacion (138 MB, ~2 min)...\n")
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250589/suppl/GSE250589_EAE_multiomics_annotation.csv.gz"
  dir.create(dirname(ann_file), showWarnings = FALSE, recursive = TRUE)
  download.file(url, ann_file, mode = "wb", quiet = FALSE)
}
cat("Archivo:", ann_file,
    " (", round(file.size(ann_file)/1024/1024, 1), " MB)\n\n")

# --- 2. Cargar con data.table (mucho mas rapido que read.csv) ---
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Instala data.table primero: install.packages('data.table')")
}
suppressPackageStartupMessages(library(data.table))

cat("Leyendo CSV comprimido...\n")
ann <- fread(ann_file, data.table = FALSE)

cat("Dimensiones:", nrow(ann), "x", ncol(ann), "\n\n")
cat("Columnas disponibles:\n")
print(colnames(ann))

cat("\nPrimeras 3 filas:\n")
print(head(ann, 3))

# --- 3. Detectar columnas clave ---
find_col <- function(patterns, cols) {
  for (p in patterns) {
    hits <- grep(p, cols, ignore.case = TRUE, value = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  return(NA)
}

col_celltype <- find_col(c("^celltype$", "cell_type", "annotation", "cluster",
                            "^class$", "cellclass", "cell_class"), colnames(ann))
col_condition <- find_col(c("^condition$", "^group$", "^stage$", "disease",
                             "^time$", "treatment", "^status$"), colnames(ann))
col_sample <- find_col(c("^sample$", "replicate", "orig.ident", "^batch$",
                          "library", "^donor$", "^mouse$", "^subject$"), colnames(ann))

cat("\n--- Columnas clave detectadas ---\n")
cat("celltype:  ", ifelse(is.na(col_celltype), "<no detectada>", col_celltype), "\n")
cat("condition: ", ifelse(is.na(col_condition), "<no detectada>", col_condition), "\n")
cat("sample:    ", ifelse(is.na(col_sample), "<no detectada>", col_sample), "\n")

# --- 4. Resumen de la estructura ---
cat("\n", SEP, "\n")
cat("ESTRUCTURA DEL DATASET\n")
cat(SEP, "\n\n")

if (!is.na(col_celltype)) {
  cat("--- Tipos celulares ---\n")
  tc <- table(ann[[col_celltype]], useNA = "ifany")
  print(tc)
  cat("Total celulas:", sum(tc), "\n\n")
}

if (!is.na(col_condition)) {
  cat("--- Condiciones ---\n")
  cond <- table(ann[[col_condition]], useNA = "ifany")
  print(cond)
  cat("\n")
}

if (!is.na(col_sample)) {
  cat("--- Replicas / muestras ---\n")
  cat("Numero de muestras unicas:", length(unique(ann[[col_sample]])), "\n")
  print(table(ann[[col_sample]], useNA = "ifany"))
  cat("\n")
}

if (!is.na(col_celltype) && !is.na(col_condition)) {
  cat("--- Crosstab celltype x condition ---\n")
  ct <- table(ann[[col_celltype]], ann[[col_condition]], useNA = "ifany")
  print(ct)
  cat("\n")
}

# --- 5. Foco en oligodendrocitos (MOL1, MOL2, MOL5/6, OPC, COP) ---
if (!is.na(col_celltype)) {
  olg_patterns <- c("MOL", "OPC", "COP", "Oligo")
  olg_vals <- unique(ann[[col_celltype]])
  olg_hits <- olg_vals[grepl(paste(olg_patterns, collapse = "|"),
                              olg_vals, ignore.case = TRUE)]
  if (length(olg_hits) > 0) {
    cat("\n", SEP, "\n")
    cat("OLIGODENDROCITOS Y PROGENITORES\n")
    cat(SEP, "\n\n")
    cat("Subtipos detectados:", paste(olg_hits, collapse = ", "), "\n\n")

    ann_olg <- ann[ann[[col_celltype]] %in% olg_hits, ]
    cat("Celulas de la linea OLG:", nrow(ann_olg), "\n\n")

    if (!is.na(col_condition)) {
      cat("--- Distribucion por condicion (solo linea OLG) ---\n")
      ct_olg <- table(ann_olg[[col_celltype]], ann_olg[[col_condition]],
                       useNA = "ifany")
      print(ct_olg)
      cat("\n")
    }

    if (!is.na(col_sample) && !is.na(col_condition)) {
      cat("--- Replicas por condicion (solo linea OLG) ---\n")
      rep_by_cond <- ann_olg %>%
        group_by(condition = .data[[col_condition]]) %>%
        summarise(n_samples = n_distinct(.data[[col_sample]]),
                  n_cells   = n(),
                  .groups   = "drop")
      print(as.data.frame(rep_by_cond))
      cat("\n")
    }

    # Guardar tabla resumen de OLG
    resumen_olg <- ann_olg %>%
      group_by(celltype = .data[[col_celltype]],
               condition = if (!is.na(col_condition)) .data[[col_condition]] else "NA") %>%
      summarise(n_cells = n(), .groups = "drop")
    write.csv(resumen_olg,
              "results/29_GSE250589_OLG_structure.csv",
              row.names = FALSE)
    cat("Guardado: results/29_GSE250589_OLG_structure.csv\n")
  }
}

# --- 6. Verificar si hay columnas de QC o metadata adicionales ---
qc_cols <- grep("nCount|nFeature|percent|pct|mito|ribo|QC|score|UMAP|PC_",
                 colnames(ann), ignore.case = TRUE, value = TRUE)
if (length(qc_cols) > 0) {
  cat("\n--- Columnas de QC / embebidos detectadas ---\n")
  cat(paste(qc_cols, collapse = ", "), "\n")
  cat("Esto sugiere que la anotacion contiene metadatos ricos (QC, UMAP).\n")
} else {
  cat("\nSin columnas de QC ni UMAP en la anotacion (solo metadata basica).\n")
}

# --- 7. Resumen para el anteproyecto ---
cat("\n", SEP, "\n")
cat("RESUMEN PARA EL ANTEPROYECTO\n")
cat(SEP, "\n\n")
cat("Dataset: GSE250589 (Zheng et al. 2023, preprint bioRxiv)\n")
cat("Modelo: EAE en raton Sox10:Cre-RCE:LoxP (linaje OLG etiquetado)\n")
cat("Tecnologia: single-cell multiome (RNA + ATAC)\n")
if (!is.na(col_celltype) && !is.na(col_condition)) {
  ct <- table(ann[[col_celltype]], ann[[col_condition]])
  olg_rows <- grep("MOL|OPC|COP", rownames(ct), ignore.case = TRUE)
  if (length(olg_rows) > 0) {
    cat("\nCelulas de la linea OLG por condicion:\n")
    print(ct[olg_rows, , drop = FALSE])
  }
}
cat("\nNOTA: la matriz de expresion (44.7 GB) NO se descarga en este\n")
cat("analisis por limitaciones de hardware (30 GB libres, 15 GB RAM).\n")
cat("El dataset se cita en el anteproyecto como 'future work' con\n")
cat("analisis pendiente para Aim 2.\n")

cat("\n", SEP, "\nFIN\n", SEP, "\n")
