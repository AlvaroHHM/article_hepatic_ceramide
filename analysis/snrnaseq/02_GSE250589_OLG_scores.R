# ============================================================
# 29b_GSE250589_OLG_scores.R
# Fase 1b: pseudobulk de scores precomputados en GSE250589
# NO requiere la matriz de expresion (44.7 GB)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

dir.create("results", showWarnings = FALSE)
sink("logs/29b_GSE250589_scores_log.txt", split = TRUE)
on.exit(sink())

cat("=== GSE250589 Fase 1b: pseudobulk de scores ===\n\n")

# --- 1. Leer anotacion ---
ann_file <- "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
ann <- fread(ann_file, data.table = FALSE)
cat("Dimensiones:", nrow(ann), "x", ncol(ann), "\n")

# --- 2. Parsear stage y condition desde orig.ident ---
ann <- ann %>%
  mutate(
    stage = case_when(
      grepl("^Early", orig.ident) ~ "Early",
      grepl("^Peak",  orig.ident) ~ "Peak",
      grepl("^Late",  orig.ident) ~ "Late",
      TRUE ~ NA_character_
    ),
    condition = case_when(
      grepl("Ctrl", orig.ident) ~ "Ctrl",
      grepl("^Early|^Peak|^Late", orig.ident) &
        !grepl("Ctrl", orig.ident) ~ "EAE",
      TRUE ~ NA_character_
    )
  )

cat("Distribucion stage x condition (todas las celulas):\n")
print(table(ann$stage, ann$condition, useNA = "ifany"))

# --- 3. Filtrar a subtipos OLG ---
olg_types <- c("MOL1", "MOL2", "MOL56", "OPC", "COP")
ann_olg <- ann %>% filter(cellType %in% olg_types)
cat("\nCelulas OLG:", nrow(ann_olg), "\n")
cat("Distribucion por subtipo x stage:\n")
print(table(ann_olg$cellType, ann_olg$stage))

# --- 4. Scores a agregar ---
score_cols <- c("ImmuneScore", "DamageScore", "InterferonScore",
                "ImmuneDamage", "ImmuneInterferon", "DamageInterferon")

# Verificar cuales existen
score_cols <- intersect(score_cols, colnames(ann_olg))
cat("\nScores presentes:", paste(score_cols, collapse = ", "), "\n")

# --- 5. Pseudobulk: media de score por sample x cellType ---
pseudo <- ann_olg %>%
  group_by(orig.ident, stage, condition, cellType) %>%
  summarise(
    n_cells = n(),
    across(all_of(score_cols), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 20)   # minimo 20 celulas por replica

cat("\nReplicas validas (n_cells >= 20):", nrow(pseudo), "\n")

# --- 6. Resumen por condition x cellType ---
resumen <- pseudo %>%
  group_by(stage, condition, cellType) %>%
  summarise(
    n_rep    = n(),
    n_cells  = sum(n_cells),
    across(all_of(score_cols), ~mean(.x, na.rm = TRUE)),
    .groups  = "drop"
  ) %>%
  arrange(stage, cellType, condition)

cat("\n", rep("=", 60), "\n")
cat("RESUMEN POR STAGE x CONDITION x CELLTYPE\n")
cat(rep("=", 60), "\n\n")
print(as.data.frame(resumen), digits = 3)

# --- 7. Delta EAE - Ctrl por stage x cellType (MOL2 y MOL56) ---
delta <- resumen %>%
  pivot_longer(all_of(score_cols), names_to = "score", values_to = "value") %>%
  pivot_wider(names_from = condition, values_from = value) %>%
  mutate(delta_EAE_vs_Ctrl = EAE - Ctrl) %>%
  filter(!is.na(delta_EAE_vs_Ctrl), is.finite(delta_EAE_vs_Ctrl), cellType %in% c("MOL2", "MOL56"))

cat("\n", rep("=", 60), "\n")
cat("DELTA EAE - CTRL POR STAGE x CELLTYPE (MOL2 y MOL56)\n")
cat(rep("=", 60), "\n\n")
print(as.data.frame(delta %>% select(stage, cellType, score, Ctrl, EAE, delta_EAE_vs_Ctrl)),
      digits = 3)

# --- 8. Falcao predictions (proporcion de estados por condicion) ---
falcao_cols <- grep("^cellTypeFalcao2018_prediction.score\\.",
                    colnames(ann_olg), value = TRUE)
if (length(falcao_cols) > 0) {
  # Estado asignado por Falcao (ya esta en cellTypeFalcao2018)
  if ("cellTypeFalcao2018" %in% colnames(ann_olg)) {
    falcao_table <- ann_olg %>%
      group_by(stage, condition, cellTypeFalcao2018) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(stage, condition) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup()

    cat("\n", rep("=", 60), "\n")
    cat("ESTADOS FALCAO2018 (proporcion por condition x stage)\n")
    cat(rep("=", 60), "\n\n")
    print(as.data.frame(falcao_table %>%
                          filter(grepl("MOL2|MOL5_6|MOL3", cellTypeFalcao2018)) %>%
                          select(stage, condition, cellTypeFalcao2018, prop) %>%
                          pivot_wider(names_from = condition, values_from = prop)),
          digits = 3)
  }
}

# --- 9. Guardar ---
write.csv(resumen, "results/29b_GSE250589_OLG_scores_summary.csv", row.names = FALSE)
write.csv(delta,   "results/29b_GSE250589_OLG_scores_delta.csv",   row.names = FALSE)
if (exists("falcao_table")) {
  write.csv(falcao_table, "results/29b_GSE250589_OLG_falcao.csv", row.names = FALSE)
}

cat("\n", rep("=", 60), "\n")
cat("Guardado:\n")
cat("  results/29b_GSE250589_OLG_scores_summary.csv\n")
cat("  results/29b_GSE250589_OLG_scores_delta.csv\n")
cat("  results/29b_GSE250589_OLG_falcao.csv\n")
cat(rep("=", 60), "\n")
