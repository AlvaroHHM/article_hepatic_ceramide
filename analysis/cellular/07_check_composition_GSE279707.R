# ============================================================
# 28_check_composition_GSE279707.R
# Verifica cambio de composición celular entre EAE y Healthy
# ============================================================

suppressPackageStartupMessages({library(dplyr)})

# Cargar los dos CSVs del script 26
pos <- read.csv("results/26_GSE279707_pos.csv")
neg <- read.csv("results/26_GSE279707_neg.csv")

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

cat("\n--- Diagnóstico de composición ---\n")
cat("Interpretación:\n")
cat("  - Si Snap25/Syp/Map2 caen en POS => pérdida de contenido neuronal\n")
cat("  - Si Gfap cae en NEG => pérdida de contenido glial\n")
cat("  - Si Irf7 sube en ambos => infiltración inmune\n")
cat("  - Cualquier combinación => log2FC confundido por composición\n")
