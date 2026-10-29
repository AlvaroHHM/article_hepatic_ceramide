# ============================================================
# 35a_diagnose_h5_and_barcodes.R (v2)
# Fix: h5$ls() sin 'depth' (hdf5r 1.3.16 no lo acepta)
# ============================================================

suppressPackageStartupMessages({
  library(hdf5r)
  library(dplyr)
})

dir.create("results", showWarnings = FALSE)
sink("results/35a_h5_structure.txt", split = TRUE)
on.exit(sink())

cat("=== H5 + BARCODE DIAGNOSTIC (v2) ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# ------------------------------------------------------------
# 1. Cargar mapa
# ------------------------------------------------------------
map <- read.csv("results/34_GSM_mapping_v2.csv", stringsAsFactors = FALSE)
map_use <- map %>% filter(!is.na(orig.ident), !grepl("Naive", timepoint))
cat("GSMs utilizables:", nrow(map_use), "\n")

# ------------------------------------------------------------
# 2. Estructura de un .h5
# ------------------------------------------------------------
h5_dir <- "data/GSE250589/h5"
first_gsm <- map_use$GSM[1]
first_file <- list.files(h5_dir, pattern = paste0("^", first_gsm, "_"),
                         full.names = TRUE)
first_file <- first_file[grepl("\\.h5$", first_file)][1]
cat("\n--- Archivo de prueba:", basename(first_file), "---\n")

h5 <- H5File$new(first_file, mode = "r")
cat("\nRutas HDF5 nivel superior:\n")
print(h5$ls()$name)
cat("\nEstructura recursiva:\n")
print(h5$ls(recursive = TRUE)$name)
h5$close_all()

# ------------------------------------------------------------
# 3. Detectar paths de features/barcodes
# ------------------------------------------------------------
cat("\n--- Formato de genes y barcodes (.h5) ---\n")
h5 <- H5File$new(first_file, mode = "r")

# Probar cellranger-arc (con/sin anidacion)
possible_feat <- c(
  "/matrix/features/name",
  "/matrix/features/id",
  "/matrix/Gene Expression/features/name",
  "/matrix/Gene Expression/features/id"
)
possible_bc <- c(
  "/matrix/barcodes",
  "/matrix/Gene Expression/barcodes"
)

features_path <- NULL; barcodes_path <- NULL
for (p in possible_feat) if (h5$exists(p)) { features_path <- p; break }
for (p in possible_bc)   if (h5$exists(p)) { barcodes_path <- p; break }

# Fallback: buscar con ls recursivo
if (is.null(features_path) || is.null(barcodes_path)) {
  paths <- h5$ls(recursive = TRUE)$name
  if (is.null(features_path)) {
    hit <- paths[grepl("features/name$", paths)]
    if (length(hit)) features_path <- hit[1]
  }
  if (is.null(barcodes_path)) {
    hit <- paths[grepl("barcodes$", paths)]
    if (length(hit)) barcodes_path <- hit[1]
  }
}

if (!is.null(features_path)) {
  features <- h5[[features_path]]$read()
  cat("Features path:", features_path, "\n")
  cat("Num features:", length(features), "\n")
  cat("Primeros 10:", paste(head(features, 10), collapse = ", "), "\n")
  is_ens <- grepl("^ENSMUSG|^ENSG|^ENS", features[1:min(20, length(features))])
  cat("Formato:", if (sum(is_ens) >= 5) "ENSEMBL" else "SYMBOL", "\n")
} else {
  cat("[WARN] No se encontro ruta de features\n")
}

if (!is.null(barcodes_path)) {
  barcodes_h5 <- h5[[barcodes_path]]$read()
  cat("\nBarcodes path:", barcodes_path, "\n")
  cat("Num barcodes:", length(barcodes_h5), "\n")
  cat("Primeros 5:", paste(head(barcodes_h5, 5), collapse = ", "), "\n")
  suffix_h5 <- unique(sub(".*(-\\d+)$", "\\1", barcodes_h5))
  cat("Sufijos en .h5:", paste(head(suffix_h5, 5), collapse = ", "), "\n")
} else {
  cat("[WARN] No se encontro ruta de barcodes\n")
}
h5$close_all()

# ------------------------------------------------------------
# 4. Barcodes del CSV
# ------------------------------------------------------------
cat("\n--- Formato de barcodes (CSV) ---\n")
ann_file <- "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
ann_head <- read.csv(gzfile(ann_file), colClasses = "character",
                     check.names = FALSE)
if (colnames(ann_head)[1] == "") colnames(ann_head)[1] <- "barcode"
barcodes_csv <- ann_head$barcode
cat("Primeros 5 CSV:", paste(head(barcodes_csv, 5), collapse = ", "), "\n")
cat("Num barcodes CSV:", length(barcodes_csv), "\n")
suffix_csv <- unique(sub(".*(-\\d+)$", "\\1", head(barcodes_csv, 1000)))
cat("Sufijos en CSV (1000 primeros):", paste(head(suffix_csv, 5), collapse = ", "), "\n")

# ------------------------------------------------------------
# 5. MOL3_EAE por muestra
# ------------------------------------------------------------
cat("\n--- MOL3_EAE por muestra (CSV) ---\n")
mol3_by_oid <- ann_head %>%
  filter(cellTypeFalcao2018 == "MOL3_EAE") %>%
  count(orig.ident, sort = TRUE)
cat("Total MOL3_EAE:", sum(mol3_by_oid$n), "\n")
print(mol3_by_oid, right = FALSE)
write.csv(mol3_by_oid, "results/35a_mol3_counts.csv", row.names = FALSE)

# ------------------------------------------------------------
# 6. Test de alineacion
# ------------------------------------------------------------
cat("\n--- Test de alineacion ---\n")
h5 <- H5File$new(first_file, mode = "r")
barcodes_h5 <- h5[[barcodes_path]]$read()
h5$close_all()

first_oid <- map_use$orig.ident[1]
bc_mol3_csv <- ann_head$barcode[
  ann_head$orig.ident == first_oid & ann_head$cellTypeFalcao2018 == "MOL3_EAE"
]

cat("GSM:", first_gsm, "| orig.ident:", first_oid, "\n")
cat("Barcodes en .h5:", length(barcodes_h5), "\n")
cat("Barcodes MOL3_EAE en CSV:", length(bc_mol3_csv), "\n\n")

# Test 1: directo
cat("Test 1 (directo):", length(intersect(barcodes_h5, bc_mol3_csv)), "\n")
# Test 2: strip sufijo
cat("Test 2 (sin sufijo):",
    length(intersect(sub("-\\d+$", "", barcodes_h5),
                     sub("-\\d+$", "", bc_mol3_csv))), "\n")
# Test 3: probar sufijos del CSV
cat("Test 3 (con sufijos del CSV):\n")
for (suf in head(suffix_csv, 5)) {
  n <- length(intersect(barcodes_h5, paste0(bc_mol3_csv, suf)))
  cat("  ", suf, ":", n, "\n")
}
# Test 4: probar sufijos en h5
cat("Test 4 (con sufijos del h5):\n")
for (suf in head(suffix_h5, 3)) {
  n <- length(intersect(paste0(sub("-\\d+$","",bc_mol3_csv), suf), barcodes_h5))
  cat("  ", suf, ":", n, "\n")
}

# ------------------------------------------------------------
# 7. Estimacion de recursos
# ------------------------------------------------------------
cat("\n--- Estimacion de recursos ---\n")
sizes_mb <- sapply(map_use$GSM, function(g) {
  f <- list.files(h5_dir, pattern = paste0("^", g, "_"), full.names = TRUE)
  f <- f[grepl("\\.h5$", f)][1]
  if (length(f)) file.info(f)$size / 1e6 else NA
})
cat("Tamano .h5: min=", round(min(sizes_mb, na.rm=TRUE), 1), "MB",
    " max=", round(max(sizes_mb, na.rm=TRUE), 1), "MB",
    " mediana=", round(median(sizes_mb, na.rm=TRUE), 1), "MB\n")
cat("RAM total (cargar 26 a la vez): ~",
    round(sum(sizes_mb, na.rm=TRUE) / 1024, 1), "GB\n")

cat("\n=== FIN ===\n")
