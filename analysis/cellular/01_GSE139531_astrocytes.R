# ============================================================
# 23_reanalysis_GSE139531.R
# Reanalisis de GSE139531 (Chao et al. 2019, Cell)
# Datos de esfingolipidos en astrocitos purificados de EAE
# Objetivo: cerrar el nodo 7 (SNC) con datos celula-especificos
# Uso: Rscript scripts/23_reanalysis_GSE139531.R
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})

dir.create("data/GSE139531", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/23_GSE139531_log.txt", split = TRUE)
on.exit(sink())

message("=== Reanalisis GSE139531: esfingolipidos en astrocitos ===")
message("Fecha: ", Sys.time())

# ------------------------------------------------------------
# 1. Localizar el archivo Excel (varias rutas posibles)
# ------------------------------------------------------------
candidatos <- c(
  "data/GSE139531/GSE139531/GSE139531_sphingolipid-data-collection.xlsx",
  "data/GSE139531/GSE139531_sphingolipid-data-collection.xlsx",
  "data/GSE139531/GSE139531_sphingolipid-data-collection.csv"
)
xlsx_file <- candidatos[file.exists(candidatos)][1]

if (is.na(xlsx_file)) {
  message("No se encontro el archivo. Buscando en data/GSE139531/...")
  all_files <- list.files("data/GSE139531", pattern = "sphingolipid",
                           recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) stop("No se encontro el archivo sphingolipid-data-collection")
  xlsx_file <- all_files[1]
}

message("Archivo: ", xlsx_file)
message("Tamano: ", round(file.size(xlsx_file) / 1024 / 1024, 1), " MB")

# ------------------------------------------------------------
# 2. Listar hojas del Excel
# ------------------------------------------------------------
sheets <- excel_sheets(xlsx_file)
message("\nHojas encontradas: ", length(sheets))
for (i in seq_along(sheets)) message("  ", i, ". ", sheets[i])

# ------------------------------------------------------------
# 3. Explorar cada hoja (diagnostico)
# ------------------------------------------------------------
explorar_hoja <- function(sheet_name, n_max = 5) {
  message("\n========================================")
  message("Hoja: ", sheet_name)
  message("========================================")

  df <- tryCatch(
    read_excel(xlsx_file, sheet = sheet_name, n_max = n_max),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(df) || nrow(df) == 0) {
    message("  Hoja vacia o no leible")
    return(NULL)
  }

  message("Dimensiones (primeras ", nrow(df), " filas): ",
          paste(dim(df), collapse = " x "))
  message("\nColumnas:")
  print(colnames(df))

  message("\nPrimeras filas:")
  print(as.data.frame(df))

  # Detectar genes de interes en cualquier columna
  genes_interes <- c("B4galt6", "B4galt5", "Ugcg", "Smpd1", "Smpd3",
                     "Sgms1", "Sgms2", "Cers2", "Cers5", "Cers6",
                     "Asah1", "Sgpl1", "S1pr3", "Fa2h", "Pla2g4a", "Mavs")

  # Revisar cada columna character por genes
  for (col in colnames(df)) {
    if (is.character(df[[col]])) {
      hits <- intersect(genes_interes, unique(df[[col]]))
      if (length(hits) > 0) {
        message("\n  Columna '", col, "' contiene genes de interes: ",
                paste(hits, collapse = ", "))
      }
    }
  }

  return(df)
}

for (s in sheets) {
  explorar_hoja(s)
}

# ------------------------------------------------------------
# 4. Busqueda automatica de la hoja con datos de expresion
# ------------------------------------------------------------
# Criterios: hoja con genes de interes Y grupos (Control, EAE, KD)
genes_interes <- c("B4galt6", "B4galt5", "Ugcg", "Smpd1", "Smpd3",
                   "Sgms1", "Sgms2", "Cers2", "Cers5", "Cers6",
                   "Asah1", "Sgpl1", "S1pr3", "Fa2h", "Pla2g4a", "Mavs")

mejor_hoja <- NULL
mejor_score <- 0

for (s in sheets) {
  df <- tryCatch(
    suppressMessages(read_excel(xlsx_file, sheet = s)),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) next

  # Score: cuantos genes de interes aparecen en cualquier columna character
  score <- 0
  for (col in colnames(df)) {
    if (is.character(df[[col]])) {
      score <- score + length(intersect(genes_interes, unique(df[[col]])))
    }
  }

  message("  Hoja '", s, "': score = ", score, " genes de interes")

  if (score > mejor_score) {
    mejor_score <- score
    mejor_hoja <- s
  }
}

if (is.null(mejor_hoja) || mejor_score == 0) {
  message("\nADVERTENCIA: no se encontro hoja con genes de interes.")
  message("Revisar el log y adaptar manualmente.")
} else {
  message("\n*** Hoja seleccionada: '", mejor_hoja, "' (score = ",
          mejor_score, ") ***")
}

# ------------------------------------------------------------
# 5. Cargar la hoja completa y extraer datos
# ------------------------------------------------------------
if (!is.null(mejor_hoja)) {
  message("\n========================================")
  message("ANALISIS DE LA HOJA: ", mejor_hoja)
  message("========================================")

  df_full <- read_excel(xlsx_file, sheet = mejor_hoja)
  message("Dimensiones completas: ", paste(dim(df_full), collapse = " x "))

  # Guardar la hoja completa como CSV para inspeccion manual
  write.csv(df_full,
            paste0("results/23_GSE139531_", gsub("[^A-Za-z0-9]", "_", mejor_hoja), ".csv"),
            row.names = FALSE)
  message("Hoja completa guardada en results/")

  # Intentar encontrar columnas de genes, muestras y grupos
  # (esto depende de la estructura del Excel; puede requerir ajuste manual)

  # Estrategia 1: hoja en formato wide (genes en filas, muestras en columnas)
  col_genes <- colnames(df_full)[1]
  message("\nAsumiendo columna de genes: '", col_genes, "'")

  genes_presentes <- intersect(genes_interes, as.character(df_full[[col_genes]]))
  message("Genes de interes presentes: ", length(genes_presentes), "/",
          length(genes_interes))
  message("Genes presentes: ", paste(genes_presentes, collapse = ", "))
  message("Genes ausentes: ", paste(setdiff(genes_interes, genes_presentes),
                                     collapse = ", "))

  if (length(genes_presentes) > 0) {
    # Filtrar filas de genes de interes
    df_genes <- df_full %>%
      filter(.data[[col_genes]] %in% genes_presentes) %>%
      select(all_of(col_genes), where(is.numeric))

    write.csv(df_genes, "results/23_GSE139531_sphingolipid_genes.csv",
              row.names = FALSE)
    message("OK results/23_GSE139531_sphingolipid_genes.csv")
  }

  # Imprimir resumen de valores
  message("\nResumen de la hoja:")
  print(summary(df_full))
}

# ------------------------------------------------------------
# 6. Buscar el RNA-seq crudo (por si el Excel no tiene counts)
# ------------------------------------------------------------
message("\n========================================")
message("NOTA SOBRE DATOS CRUDOS")
message("========================================")
message("El archivo 'sphingolipid-data-collection.xlsx' es un dataset curado de esfingolipidos.")
message("Si se requiere el RNA-seq crudo de GSE139531, buscar en SRA:")
message("  - Accesion: SRP226905 (SRA study de GSE139531)")
message("  - O usar recount3: disponible_projects(organism = 'mouse') y buscar GSE139531")
message("  - O descargar de GEO: GSE139531_RAW.tar (puede contener counts)")

# Intentar acceso via recount3
message("\nIntentando acceso via recount3...")
tryCatch({
  library(recount3)
  proyectos_mouse <- available_projects(organism = "mouse")
  hit <- proyectos_mouse[grepl("139531|Chao|astrocyte|cPLA2", proyectos_mouse$project,
                                ignore.case = TRUE), ]
  if (nrow(hit) > 0) {
    message("Proyectos encontrados en recount3:")
    print(hit[, c("project", "study_title", "n_samples")])
  } else {
    message("GSE139531 no esta en recount3. Usar SRA (SRP226905) o GEO RAW.tar.")
  }
}, error = function(e) {
  message("recount3 no disponible: ", conditionMessage(e))
})

message("\n=== FIN ===")
message("Revisar el log completo: logs/23_GSE139531_log.txt")
