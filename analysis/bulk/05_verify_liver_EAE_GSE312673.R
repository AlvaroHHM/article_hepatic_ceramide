# ============================================================
# 33_verify_liver_EAE_GSE312673.R
# Analisis del transcriptoma hepatico en EAE (GSE312673)
# Contrastes: peak vs control, pre-onset vs control
# Panel: SMPD3, SGMS1, SGMS2, B4GALT5/6, CERS2, CYP7A1, NR1H4, etc.
#
# Uso: Rscript scripts/33_verify_liver_EAE_GSE312673.R
# Salidas:
#   results/33_liver_EAE_log2FC.csv
#   results/33_liver_EAE_panel.csv
#   logs/33_verify_liver_EAE_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(GEOquery)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

dir.create("data/GSE312673", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

sink("logs/33_verify_liver_EAE_log.txt", split = TRUE)
on.exit(sink())

cat("============================================================\n")
cat("GSE312673: LIVER TRANSCRIPTOME IN EAE\n")
cat("Date:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# ============================================================
# 1. DESCARGA DE ARCHIVOS SUPLEMENTARIOS
# ============================================================
cat("--- 1. Descarga de archivos suplementarios ---\n")

# getGEOSuppFiles descarga a data/GSE312673/
# El archivo suele ser un .txt.gz o .csv.gz con counts
supp_files <- tryCatch({
  getGEOSuppFiles("GSE312673", baseDir = "data/GSE312673",
                  makeDirectory = FALSE, fetch_files = TRUE)
}, error = function(e) {
  cat("  Error en getGEOSuppFiles:", e$message, "\n")
  NULL
})

if (is.null(supp_files) || nrow(supp_files) == 0) {
  cat("  No se pudieron descargar archivos suplementarios.\n")
  cat("  Alternativa manual:\n")
  cat("    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312673/suppl/\n")
  stop("Descarga manual requerida.")
}

cat("  Archivos descargados:\n")
print(supp_files)

# ============================================================
# 2. LECTURA DE LA MATRIZ DE COUNTS
# ============================================================
cat("\n--- 2. Lectura de la matriz de counts ---\n")

# Buscar el archivo de counts (no la series matrix)
count_files <- list.files("data/GSE312673", full.names = TRUE,
                          pattern = "\\.(txt|tsv|csv)(\\.gz)?$")
count_files <- count_files[!grepl("series_matrix", count_files)]

if (length(count_files) == 0) stop("No se encontro archivo de counts.")

count_file <- count_files[1]
cat("  Usando:", count_file, "\n")

# Intentar leer con read.delim (auto-detecta separador)
counts <- tryCatch({
  read.delim(count_file, row.names = 1, check.names = FALSE,
             stringsAsFactors = FALSE)
}, error = function(e) {
  # Fallback a read.csv
  read.csv(count_file, row.names = 1, check.names = FALSE,
           stringsAsFactors = FALSE)
})

cat("  Dimensiones:", paste(dim(counts), collapse = " x "), "\n")
cat("  Primeras filas/columnas:\n")
print(counts[1:3, 1:min(5, ncol(counts))])

# ============================================================
# 3. DETECCION DE COLUMNAS Y GRUPOS
# ============================================================
cat("\n--- 3. Deteccion de columnas y grupos ---\n")

# Los nombres de columna suelen incluir "peak", "pre", "control", "EAE"
samples <- colnames(counts)
cat("  Muestras:", paste(samples, collapse = ", "), "\n")

# Clasificar cada muestra en uno de 4 grupos
group <- case_when(
  grepl("peak.*EAE|EAE.*peak", samples, ignore.case = TRUE) ~ "Peak_EAE",
  grepl("peak.*(NC|control|Ctrl)|(NC|control|Ctrl).*peak",
        samples, ignore.case = TRUE) ~ "Peak_Ctrl",
  grepl("pre.*EAE|EAE.*pre", samples, ignore.case = TRUE) ~ "Pre_EAE",
  grepl("pre.*(NC|control|Ctrl)|(NC|control|Ctrl).*pre",
        samples, ignore.case = TRUE) ~ "Pre_Ctrl",
  TRUE ~ NA_character_
)

cat("  Asignacion de grupos:\n")
print(table(group, useNA = "ifany"))

if (sum(is.na(group)) > 0) {
  cat("  WARN: algunas muestras no se asignaron. Revisar nombres.\n")
}

# ============================================================
# 4. MAPEO ENSEMBL -> SYMBOL
# ============================================================
cat("\n--- 4. Mapeo Ensembl -> Symbol ---\n")

ids <- rownames(counts)
# Strip de version (.N)
ids_clean <- sub("\\..*$", "", ids)

# Detectar si son Ensembl
is_ens <- grepl("^ENSMUSG", ids_clean)
cat("  IDs Ensembl:", sum(is_ens), "/", length(ids_clean), "\n")

if (sum(is_ens) > 0) {
  symbols <- tryCatch({
    suppressMessages(
      mapIds(org.Mm.eg.db::org.Mm.eg.db,
             keys = ids_clean[is_ens],
             keytype = "ENSEMBL",
             column = "SYMBOL",
             multiVals = "first")
    )
  }, error = function(e) NULL)

  if (!is.null(symbols)) {
    # Crear vector de simbolos (NA para no-Ensembl)
    sym_vec <- rep(NA_character_, length(ids_clean))
    sym_vec[is_ens] <- as.character(symbols)
    counts$symbol <- sym_vec
    cat("  Mapeados:", sum(!is.na(sym_vec)), "/", length(sym_vec), "\n")
  } else {
    counts$symbol <- ids_clean
    cat("  WARN: mapIds fallo. Usando IDs crudos.\n")
  }
} else {
  counts$symbol <- ids_clean
  cat("  IDs no son Ensembl. Usando IDs crudos como simbolo.\n")
}

# Filtrar genes con simbolo
counts_sym <- counts[!is.na(counts$symbol) & nzchar(counts$symbol), ]
# Agregar duplicados por simbolo (suma de counts)
counts_agg <- counts_sym %>%
  group_by(symbol) %>%
  summarise(across(all_of(samples), ~sum(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  column_to_rownames("symbol")

cat("  Genes unicos tras agregacion:", nrow(counts_agg), "\n")

# ============================================================
# 5. CALCULO DE log2FC (CONTRASTE SIMPLE)
# ============================================================
cat("\n--- 5. Calculo de log2FC ---\n")

# Funcion: media de CPM por grupo
calc_cpm <- function(mat) {
  t(t(mat) / colSums(mat) * 1e6)
}

cpm <- calc_cpm(counts_agg)

# Medias por grupo
grp_levels <- unique(na.omit(group))
if (length(grp_levels) == 0) stop("No hay grupos asignados.")

mean_by_group <- sapply(grp_levels, function(g) {
  cols <- which(group == g)
  if (length(cols) == 0) return(rep(NA, nrow(cpm)))
  rowMeans(cpm[, cols, drop = FALSE], na.rm = TRUE)
})

# Contrastes: Peak_EAE vs Peak_Ctrl, Pre_EAE vs Pre_Ctrl
contrasts <- list(
  Peak_EAE_vs_Ctrl = c("Peak_EAE", "Peak_Ctrl"),
  Pre_EAE_vs_Ctrl  = c("Pre_EAE",  "Pre_Ctrl")
)

log2fc_df <- data.frame(gene = rownames(cpm))

for (nm in names(contrasts)) {
  eae_grp <- contrasts[[nm]][1]
  ctrl_grp <- contrasts[[nm]][2]

  if (!(eae_grp %in% colnames(mean_by_group)) ||
      !(ctrl_grp %in% colnames(mean_by_group))) {
    cat(sprintf("  SKIP %s: grupo faltante\n", nm))
    next
  }

  lfc <- log2(mean_by_group[, eae_grp] + 1) -
         log2(mean_by_group[, ctrl_grp] + 1)
  log2fc_df[[nm]] <- lfc
}

# ============================================================
# 6. EXTRAER PANEL DE GENES DE INTERES
# ============================================================
cat("\n--- 6. Panel de genes de interes ---\n")

panel_genes <- c(
  # Ciclo SM
  "Smpd3", "Sgms1", "Sgms2",
  # Lactosilceramida
  "B4galt5", "B4galt6", "Ugcg",
  # Sintesis de novo
  "Cers2", "Cers5", "Cers6",
  # Acidos biliares / FXR
  "Cyp7a1", "Cyp8b1", "Nr1h4", "Abcb11", "Slc10a1",
  # Inflamacion
  "Tnf", "Il6", "Il1b", "Ccl2",
  # S1P
  "Sphk1", "S1pr1"
)


cat("  DIAGNOSTICO: muestra de simbolos en log2fc_df$gene:\n")
cat("   ", paste(head(log2fc_df$gene, 10), collapse = ", "), "\n")
cat("  Buscando especificamente: Smpd3, Sgms1, B4galt6, Cyp7a1...\n")
for (g in c("Smpd3","Sgms1","B4galt6","Cyp7a1","Nr1h4")) {
  cat(sprintf("    %-8s : %s\n", g,
              if (g %in% log2fc_df$gene) "PRESENTE" else "AUSENTE"))
}

panel_present <- intersect(panel_genes, log2fc_df$gene)
cat("  Presentes:", length(panel_present), "/", length(panel_genes), "\n")
cat("  Ausentes:", paste(setdiff(panel_genes, panel_present), collapse = ", "), "\n")

panel_df <- log2fc_df %>%
  filter(gene %in% panel_present) %>%
  arrange(gene)

cat("\n  log2FC del panel:\n")
print(panel_df)

# ============================================================
# 7. GUARDAR RESULTADOS
# ============================================================
write.csv(log2fc_df, "results/33_liver_EAE_log2FC.csv", row.names = FALSE)
write.csv(panel_df, "results/33_liver_EAE_panel.csv", row.names = FALSE)

cat("\n  -> results/33_liver_EAE_log2FC.csv\n")
cat("  -> results/33_liver_EAE_panel.csv\n")

# ============================================================
# 8. INTERPRETACION RAPIDA
# ============================================================
cat("\n============================================================\n")
cat("INTERPRETACION\n")
cat("============================================================\n")

if ("Smpd3" %in% panel_df$gene) {
  smpd3 <- panel_df %>% filter(gene == "Smpd3")
  cat(sprintf("  Smpd3 en Peak EAE vs Ctrl: log2FC = %+.2f\n",
              smpd3$Peak_EAE_vs_Ctrl))
  cat(sprintf("  Smpd3 en Pre EAE vs Ctrl:  log2FC = %+.2f\n",
              smpd3$Pre_EAE_vs_Ctrl))
}

if ("Sgms1" %in% panel_df$gene) {
  sgms1 <- panel_df %>% filter(gene == "Sgms1")
  cat(sprintf("  Sgms1 en Peak EAE vs Ctrl: log2FC = %+.2f\n",
              sgms1$Peak_EAE_vs_Ctrl))
}

if ("B4galt6" %in% panel_df$gene) {
  b4galt6 <- panel_df %>% filter(gene == "B4galt6")
  cat(sprintf("  B4galt6 en Peak EAE vs Ctrl: log2FC = %+.2f\n",
              b4galt6$Peak_EAE_vs_Ctrl))
}

cat("\n  Comparar con el patron de Mdr2-KO (higado cronico):\n")
cat("    Mdr2-KO: Smpd3 +4.80 (padj 2e-58), Sgms1 +0.01 (NS)\n")
cat("    EAE (este analisis): ver valores arriba\n")
cat("\n  Si Smpd3 se induce en EAE Y Sgms1 permanece plano,\n")
cat("  el patron hepatico se confirma en el modelo de enfermedad.\n")
cat("\n=== FIN ===\n")
