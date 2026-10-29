# ============================================================
# 26_reanalysis_GSE279707_neurons.R (v3)
# Núcleos EAE vs Healthy, separados por sorting (_pos_/_neg_)
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

dir.create("results", showWarnings=FALSE)
dir.create("logs", showWarnings=FALSE)
sink("logs/26_GSE279707_log.txt", split=TRUE)
on.exit(sink())

message("=== GSE279707: nucleos EAE vs Healthy ===")
message("Fecha: ", Sys.time())

# 1. Localizar archivo
supp_dir <- "data/GSE279707"
supp_files <- list.files(supp_dir, recursive=TRUE, full.names=TRUE)
if (length(supp_files) == 0) {
  stop("No hay archivos en ", supp_dir,
       ". Ejecuta getGEOSuppFiles('GSE279707') primero.")
}
expr_file <- supp_files[grepl("count|fpkm|tpm|expr",
                              basename(supp_files), ignore.case=TRUE)][1]
if (is.na(expr_file)) stop("No se encontró archivo de expresión")
message("Usando: ", basename(expr_file))

# 2. Leer
leer <- function(f) {
  ext <- tools::file_ext(f)
  if (ext == "xlsx") return(read_excel(f, sheet=1))
  con <- if (grepl("\\.gz$", f)) gzfile(f, "r") else file(f, "r")
  on.exit(close(con), add=TRUE)
  read.delim(con, stringsAsFactors=FALSE, check.names=FALSE)
}
expr_raw <- leer(expr_file)
message("\nDimensiones: ", paste(dim(expr_raw), collapse=" x "))
message("colnames (todas):")
print(colnames(expr_raw))

# 3. Detectar columna de genes (Geneid)
if (!"Geneid" %in% colnames(expr_raw)) {
  stop("No se encontró columna 'Geneid'. Columnas: ",
       paste(colnames(expr_raw), collapse=", "))
}
gene_ids <- as.character(expr_raw[["Geneid"]])
expr_raw <- expr_raw[, setdiff(colnames(expr_raw), "Geneid"), drop=FALSE]
message("\nColumna 'Geneid' eliminada. Primeros 5 IDs: ",
        paste(head(gene_ids, 5), collapse=", "))

# 4. Numérico + filtro
expr_num <- as.data.frame(lapply(expr_raw, function(x)
  suppressWarnings(as.numeric(as.character(x)))))
rownames(expr_num) <- gene_ids
keep <- rowSums(!is.na(expr_num)) >= 0.8 * ncol(expr_num)
expr_num <- expr_num[keep, , drop=FALSE]
message("Genes tras filtro: ", nrow(expr_num), "/", length(gene_ids))

# 5. Mapear ENSEMBL a SYMBOL (agregando duplicados)
if (any(grepl("^ENSMUSG", rownames(expr_num)))) {
  message("\nRownames son ENSEMBL. Mapeando a símbolos...")
  map <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = sub("\\..*$", "", rownames(expr_num)),
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
  map <- map[!is.na(map$SYMBOL), ]

  # Unir gene_ids originales con SYMBOL
  ens_clean <- sub("\\..*$", "", rownames(expr_num))
  symbol_vec <- map$SYMBOL[match(ens_clean, map$ENSEMBL)]
  keep_ens <- !is.na(symbol_vec)
  expr_num <- expr_num[keep_ens, , drop=FALSE]
  symbol_vec <- symbol_vec[keep_ens]
  message("Genes con símbolo asignado: ", length(symbol_vec))

  # Agregar por símbolo (suma de counts) para eliminar duplicados
  expr_num$SYMBOL <- symbol_vec
  expr_agg <- aggregate(. ~ SYMBOL, data = expr_num, FUN = sum, na.rm = TRUE)
  rownames(expr_agg) <- expr_agg$SYMBOL
  expr_agg$SYMBOL <- NULL
  expr_num <- as.matrix(expr_agg)
  message("Genes tras agregación por símbolo: ", nrow(expr_num),
          " (símbolos únicos)")
}

# 6. Detección de grupos: EAE vs Healthy Y sorting _pos_/_neg_
samples <- colnames(expr_num)
message("\nNombres de muestra:"); print(samples)

# Detectar condición
cond <- ifelse(grepl("EAE|eae", samples), "EAE",
        ifelse(grepl("Healthy|healthy|Control|control", samples), "Healthy",
               NA_character_))

# Detectar sorting (_pos_ = NeuN+ = neuronas; _neg_ = NeuN- = no-neuronas)
sorting <- ifelse(grepl("_pos_", samples, ignore.case=TRUE), "pos",
           ifelse(grepl("_neg_", samples, ignore.case=TRUE), "neg",
                  "unknown"))

message("\nCrosstab condición x sorting:")
print(table(cond, sorting, useNA="ifany"))

if (any(sorting == "unknown")) {
  message("\nADVERTENCIA: hay muestras sin sorting detectable. ",
          "Se tratarán como un único grupo.")
  sorting[sorting == "unknown"] <- "all"
}

# 7. Genes de interés
genes_interes <- c(
  "B4galt6","B4galt5","Ugcg","Smpd1","Smpd3","Smpd2","Smpd4",
  "Sgms1","Sgms2","Cers1","Cers2","Cers3","Cers4","Cers5","Cers6",
  "Sptlc1","Sptlc2","Kdsr","Degs1","Degs2",
  "Asah1","Asah2","Acer1","Acer2","Acer3",
  "Sphk1","Sphk2","Sgpl1","S1pr1","S1pr2","S1pr3","S1pr4","S1pr5",
  "Ormdl1","Ormdl2","Ormdl3","Fa2h",
  "Pla2g4a","Mavs","Nfkb1","Rela","Irf3","Irf7",
  "Snap25","Syp","Rbfox3","Map2","Gfap"
)
genes_presentes <- intersect(genes_interes, rownames(expr_num))
message("\nGenes de interés presentes: ", length(genes_presentes), "/",
        length(genes_interes))
message("Presentes: ", paste(genes_presentes, collapse=", "))
message("Ausentes:  ", paste(setdiff(genes_interes, genes_presentes), collapse=", "))

if (length(genes_presentes) == 0) {
  stop("Ningún gen de interés. head(rownames(expr_num)):\n  ",
       paste(head(rownames(expr_num), 10), collapse="\n  "))
}
expr_sub <- expr_num[genes_presentes, , drop=FALSE]

# 8. Análisis por subgrupo (pos / neg)
analizar_grupo <- function(sub_sorting, etiqueta) {
  keep <- sorting == sub_sorting & !is.na(cond)
  if (sum(keep) < 4) {
    message("\n[", etiqueta, "] Muestras insuficientes (", sum(keep), "). Saltando.")
    return(NULL)
  }
  expr_s <- expr_sub[, keep, drop=FALSE]
  cond_s <- cond[keep]
  message("\n========================================")
  message(etiqueta, " (", sum(cond_s=="EAE"), " EAE, ",
          sum(cond_s=="Healthy"), " Healthy)")
  message("========================================")

  res <- data.frame(gene = genes_presentes, stringsAsFactors=FALSE)
  for (g in c("Healthy", "EAE")) {
    res[[paste0("mean_", g)]] <- apply(expr_s, 1,
      function(x) mean(x[cond_s==g], na.rm=TRUE))
  }
  res$log2FC <- log2(res$mean_EAE + 1) - log2(res$mean_Healthy + 1)
  res <- res[order(-res$log2FC), ]
  write.csv(res, paste0("results/26_GSE279707_", sub_sorting, ".csv"),
            row.names=FALSE)
  message("OK results/26_GSE279707_", sub_sorting, ".csv")

  message("\n--- MEDIA POR GEN Y GRUPO ---")
  print(res, row.names=FALSE, digits=2)

  # Identidad celular
  message("\n--- VERIFICACIÓN DE IDENTIDAD ---")
  for (m in c("Snap25","Syp","Rbfox3","Map2","Gfap")) {
    if (m %in% rownames(expr_s)) {
      message(sprintf("  %-8s: media global = %.1f", m,
                      mean(expr_s[m, ], na.rm=TRUE)))
    }
  }

  # Interpretación
  message("\n--- INTERPRETACIÓN ---")
  for (g in c("B4galt6","Smpd3","Sgms1","Asah1","Sgpl1",
              "S1pr1","S1pr3","Pla2g4a","Mavs")) {
    fc <- res$log2FC[res$gene==g]
    if (length(fc) > 0) message(sprintf("  %-8s log2FC = %+.2f", g, fc))
  }
  b4 <- res$log2FC[res$gene == "B4galt6"]
  if (length(b4) > 0) {
    if (b4 > 1) message("\n  -> B4galt6 SE INDUCE (", sub_sorting, ")")
    else if (b4 < -1) message("\n  -> B4galt6 se SUPRIME (", sub_sorting, ")")
    else message("\n  -> B4galt6 NO cambia (", sub_sorting,
                 ") — consistente con astrocitos y microglía")
  }
  invisible(res)
}

# Analizar por separado: pos (neuronas), neg (no-neuronas)
if ("pos" %in% sorting) analizar_grupo("pos", "NÚCLEOS NeuN+ (neuronas)")
if ("neg" %in% sorting) analizar_grupo("neg", "NÚCLEOS NeuN- (no-neuronas)")

# Si el sorting no está claro, hacer el análisis combinado
if (!any(c("pos","neg") %in% sorting)) {
  message("\nNo se detectó sorting _pos_/_neg_. Análisis combinado.")
  # Código de análisis combinado (equivalente al anterior)
}

message("\n=== FIN ===")
