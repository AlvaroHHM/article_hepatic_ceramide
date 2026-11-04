# ============================================================
# 07_check_composition_GSE279707.R
# Verifica cambio de composición celular entre EAE y Healthy
# ============================================================

suppressPackageStartupMessages({library(dplyr)})

dir.create("results/cellular", showWarnings = FALSE, recursive = TRUE)
dir.create("logs/cellular", showWarnings = FALSE, recursive = TRUE)
sink("logs/cellular/07_check_composition_GSE279707_log.txt", split = TRUE)
on.exit(sink())

cat("=== 07_check_composition_GSE279707.R ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# Cargar los dos CSVs del script 04 (antes 26_*)
pos <- read.csv("results/cellular/04_GSE279707_pos.csv", stringsAsFactors = FALSE)
neg <- read.csv("results/cellular/04_GSE279707_neg.csv", stringsAsFactors = FALSE)
if (nrow(pos) == 0 || nrow(neg) == 0) {
  stop("Faltan CSVs de 04. Ejecuta: Rscript analysis/cellular/04_GSE279707_neurons.R")
}

# Marcadores de identidad
neuronal <- c("Snap25", "Syp", "Rbfox3", "Map2")
glial    <- c("Gfap", "Aqp4", "Mbp", "Cspg4")  # Aqp4/Cspg4 si están
immune   <- c("Ptprc", "Cd68", "Cd3e", "Irf7")

check_markers <- function(df, name) {
  cat("\n---", name, "---\n")
  for (m in c(neuronal, glial, immune)) {
    if (m %in% df$gene) {
      row <- df[df$gene == m, ]
      cat(sprintf("  %-8s: Healthy=%8.1f  EAE=%8.1f  log2FC=%+.2f\n",
                  m, row$mean_Healthy, row$mean_EAE, row$log2FC))
    }
  }
}

check_markers(pos, "POS (neuronas enriquecidas)")
check_markers(neg, "NEG (glía + infiltrado)")

# Guardar tabla de marcadores (respaldo de las líneas 857-861 del .tex)
markers_panel <- c("Snap25","Syp","Rbfox3","Map2","Gfap","Irf7")
comp_tbl <- data.frame(
  marker = markers_panel,
  pos_Healthy = sapply(markers_panel, function(m) {
    v <- pos$mean_Healthy[pos$gene == m]; if (length(v)) v[1] else NA_real_
  }),
  pos_EAE = sapply(markers_panel, function(m) {
    v <- pos$mean_EAE[pos$gene == m]; if (length(v)) v[1] else NA_real_
  }),
  pos_log2FC = sapply(markers_panel, function(m) {
    v <- pos$log2FC[pos$gene == m]; if (length(v)) v[1] else NA_real_
  }),
  neg_Healthy = sapply(markers_panel, function(m) {
    v <- neg$mean_Healthy[neg$gene == m]; if (length(v)) v[1] else NA_real_
  }),
  neg_EAE = sapply(markers_panel, function(m) {
    v <- neg$mean_EAE[neg$gene == m]; if (length(v)) v[1] else NA_real_
  }),
  neg_log2FC = sapply(markers_panel, function(m) {
    v <- neg$log2FC[neg$gene == m]; if (length(v)) v[1] else NA_real_
  })
)
write.csv(comp_tbl, "results/cellular/07_check_composition_markers.csv", row.names = FALSE)
cat("OK results/cellular/07_check_composition_markers.csv\n")

cat("\n--- Diagnóstico de composición ---\n")
cat("Interpretación:\n")
cat("  - Si Snap25/Syp/Map2 caen en POS => pérdida de contenido neuronal\n")
cat("  - Si Gfap cae en NEG => pérdida de contenido glial\n")
cat("  - Si Irf7 sube en ambos => infiltración inmune\n")
cat("  - Cualquier combinación => log2FC confundido por composición\n")
