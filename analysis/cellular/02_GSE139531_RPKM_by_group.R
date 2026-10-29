# ============================================================
# 24_GSE139531_astrocytes_expression.R
# Extrae RPKM de genes de esfingolipidos en astrocitos EAE
# (shControl) y knockdowns (shB4GALT6, shUGCG, etc.)
# Cierra el nodo 7 con datos celula-especificos de Chao 2019
# Uso: Rscript scripts/24_GSE139531_astrocytes_expression.R
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/24_GSE139531_expression_log.txt", split = TRUE)
on.exit(sink())

message("=== GSE139531: astrocitos EAE, expresion de esfingolipidos ===")
message("Fecha: ", Sys.time())

# ------------------------------------------------------------
# 1. Localizar el Excel
# ------------------------------------------------------------
xlsx_file <- c(
  "data/GSE139531/GSE139531/GSE139531_sphingolipid-data-collection.xlsx",
  "data/GSE139531/GSE139531_sphingolipid-data-collection.xlsx"
)[file.exists(c(
  "data/GSE139531/GSE139531/GSE139531_sphingolipid-data-collection.xlsx",
  "data/GSE139531/GSE139531_sphingolipid-data-collection.xlsx"
))][1]

if (is.na(xlsx_file)) stop("No se encontro el archivo Excel")
message("Archivo: ", xlsx_file)

# ------------------------------------------------------------
# 2. Leer la hoja saltando la fila "Group"
# ------------------------------------------------------------
# La hoja tiene 19 columnas: NAME + 18 muestras (G1..D11)
# Fila 1 = "Group" + nombres de muestra
# Fila 2+ = genes + RPKM
raw <- suppressMessages(
  read_excel(xlsx_file, sheet = "CC_rpkm_genes_removed0.txt")
)

# Extraer la fila "Group" (nombres de muestra)
group_row <- raw[1, , drop = FALSE]
sample_names <- as.character(group_row[1, -1])  # G1..D11 -> shControl-1, etc.
message("\nNombres de muestra (fila Group):")
print(sample_names)

# Eliminar la fila "Group" y convertir a numérico
expr_df <- raw[-1, , drop = FALSE]
gene_ids <- as.character(expr_df[[1]])
expr_df <- expr_df[, -1, drop = FALSE]
expr_df[] <- lapply(expr_df, function(x) as.numeric(as.character(x)))

# Asignar nombres de muestra a las columnas
colnames(expr_df) <- sample_names
expr_df <- cbind(gene = gene_ids, expr_df)

message("\nDimensiones tras parseo: ", paste(dim(expr_df), collapse = " x "))

# ------------------------------------------------------------
# 3. Definir grupos desde los nombres de muestra
# ------------------------------------------------------------
group <- ifelse(grepl("shControl", sample_names, ignore.case = TRUE), "shControl",
         ifelse(grepl("shB4GALT6", sample_names, ignore.case = TRUE), "shB4GALT6",
         ifelse(grepl("shUGCG", sample_names, ignore.case = TRUE), "shUGCG",
         ifelse(grepl("shCERS5", sample_names, ignore.case = TRUE), "shCERS5",
         ifelse(grepl("shCERS6", sample_names, ignore.case = TRUE), "shCERS6",
         ifelse(grepl("shPLA2G4A", sample_names, ignore.case = TRUE), "shPLA2G4A",
                NA_character_))))))
names(group) <- sample_names

message("\nGrupos detectados:")
print(table(group))

# ------------------------------------------------------------
# 4. Genes de interes
# ------------------------------------------------------------
genes_interes <- c(
  # Via B4GALT/UGCG (LacCer)
  "B4galt6", "B4galt5", "Ugcg",
  # Hidrolisis de esfingomielina
  "Smpd1", "Smpd3", "Smpd2", "Smpd4",
  # SMS (sintesis de esfingomielina)
  "Sgms1", "Sgms2",
  # Sintesis de novo
  "Cers1", "Cers2", "Cers3", "Cers4", "Cers5", "Cers6",
  "Sptlc1", "Sptlc2", "Kdsr", "Degs1", "Degs2",
  # Catabolismo
  "Asah1", "Asah2", "Acer1", "Acer2", "Acer3",
  # Via S1P
  "Sphk1", "Sphk2", "Sgpl1", "S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5",
  # Regulacion
  "Ormdl1", "Ormdl2", "Ormdl3", "Fa2h",
  # cPLA2-MAVS
  "Pla2g4a", "Mavs", "Nfkb1", "Rela", "Irf3", "Irf7"
)

genes_presentes <- intersect(genes_interes, expr_df$gene)
genes_ausentes <- setdiff(genes_interes, expr_df$gene)
message("\nGenes de interes presentes: ", length(genes_presentes), "/",
        length(genes_interes))
message("Genes presentes: ", paste(genes_presentes, collapse = ", "))
message("Genes ausentes: ", paste(genes_ausentes, collapse = ", "))

# ------------------------------------------------------------
# 5. Extraer RPKM de genes de interes
# ------------------------------------------------------------
expr_genes <- expr_df %>%
  filter(gene %in% genes_presentes) %>%
  as.data.frame()

# Guardar RPKM por muestra
write.csv(expr_genes, "results/24_GSE139531_genes_RPKM_by_sample.csv",
          row.names = FALSE)
message("OK results/24_GSE139531_genes_RPKM_by_sample.csv")

# ------------------------------------------------------------
# 6. Calcular RPKM medio por grupo
# ------------------------------------------------------------
expr_long <- expr_genes %>%
  pivot_longer(-gene, names_to = "sample", values_to = "RPKM") %>%
  mutate(group = group[sample])

expr_means <- expr_long %>%
  group_by(gene, group) %>%
  summarise(
    n        = n(),
    RPKM_mean = mean(RPKM, na.rm = TRUE),
    RPKM_sd   = sd(RPKM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene, group)

write.csv(expr_means, "results/24_GSE139531_genes_RPKM_by_group.csv",
          row.names = FALSE)
message("OK results/24_GSE139531_genes_RPKM_by_group.csv")

# ------------------------------------------------------------
# 7. Tabla ancha: genes en filas, grupos en columnas
# ------------------------------------------------------------
expr_wide <- expr_means %>%
  select(gene, group, RPKM_mean) %>%
  pivot_wider(names_from = group, values_from = RPKM_mean)

message("\n========================================")
message("RPKM MEDIO POR GEN Y GRUPO (astrocitos EAE)")
message("========================================")
print(as.data.frame(expr_wide), row.names = FALSE, digits = 2)

# ------------------------------------------------------------
# 8. Log2FC de cada knockdown vs shControl
# ------------------------------------------------------------
message("\n========================================")
message("LOG2FC DE CADA KNOCKDOWN vs shControl (EAE)")
message("========================================")

control_means <- expr_means %>%
  filter(group == "shControl") %>%
  select(gene, ctrl = RPKM_mean)

knockdowns <- setdiff(unique(expr_means$group), "shControl")

log2fc_df <- lapply(knockdowns, function(kd) {
  kd_means <- expr_means %>%
    filter(group == kd) %>%
    select(gene, kd = RPKM_mean)
  inner_join(control_means, kd_means, by = "gene") %>%
    mutate(
      log2FC = log2(kd + 1) - log2(ctrl + 1),
      knockdown = kd
    ) %>%
    select(gene, knockdown, ctrl, kd, log2FC)
}) %>% bind_rows()

write.csv(log2fc_df, "results/24_GSE139531_log2FC_vs_shControl.csv",
          row.names = FALSE)
message("OK results/24_GSE139531_log2FC_vs_shControl.csv")

# Mostrar solo el knockdown de B4galt6
message("\n--- shB4GALT6 vs shControl (validacion del knockdown) ---")
kd_b4 <- log2fc_df %>%
  filter(knockdown == "shB4GALT6") %>%
  arrange(desc(abs(log2FC)))
print(as.data.frame(kd_b4), row.names = FALSE, digits = 2)

# ------------------------------------------------------------
# 9. Interpretacion para el nodo 7
# ------------------------------------------------------------
message("\n========================================")
message("INTERPRETACION PARA EL NODO 7")
message("========================================")

# a) Expresion basal de B4galt6 en shControl (astrocitos EAE)
b4_ctrl <- expr_means %>%
  filter(gene == "B4galt6", group == "shControl") %>%
  pull(RPKM_mean)

if (length(b4_ctrl) == 1) {
  message("\nB4galt6 en shControl (astrocitos EAE): RPKM = ",
          round(b4_ctrl, 1))
  if (b4_ctrl > 20) {
    message("  -> B4galt6 esta EXPRESADO en astrocitos EAE (RPKM > 20)")
    message("  -> Consistente con Mayo 2014 y Chao 2019: B4GALT6 se induce en astrocitos")
  } else if (b4_ctrl > 5) {
    message("  -> B4galt6 tiene expresion MODERADA en astrocitos EAE (5-20 RPKM)")
  } else {
    message("  -> B4galt6 tiene expresion BAJA en astrocitos EAE (<5 RPKM)")
    message("  -> Revisar: el bulk de SRP071304 lo vio plano, pero en astrocitos podria ser diferente")
  }
}

# b) Verificar que el knockdown funciono
b4_kd <- expr_means %>%
  filter(gene == "B4galt6", group == "shB4GALT6") %>%
  pull(RPKM_mean)

if (length(b4_kd) == 1 && length(b4_ctrl) == 1) {
  fold_reduction <- b4_ctrl / max(b4_kd, 0.01)
  message("\nB4galt6 en shB4GALT6: RPKM = ", round(b4_kd, 1))
  message("Reduccion por knockdown: ", round(fold_reduction, 1), "x")
  if (fold_reduction > 2) {
    message("  -> El knockdown de B4GALT6 funciono (reduccion > 2x)")
  } else {
    message("  -> ADVERTENCIA: el knockdown no redujo B4galt6 de forma clara")
  }
}

# c) Comparar SMPD3 vs SMS en astrocitos
smpd3_ctrl <- expr_means %>%
  filter(gene == "Smpd3", group == "shControl") %>%
  pull(RPKM_mean)
sgms1_ctrl <- expr_means %>%
  filter(gene == "Sgms1", group == "shControl") %>%
  pull(RPKM_mean)
sgms2_ctrl <- expr_means %>%
  filter(gene == "Sgms2", group == "shControl") %>%
  pull(RPKM_mean)

if (length(smpd3_ctrl) == 1) {
  message("\nEn astrocitos EAE (shControl):")
  message("  Smpd3 RPKM = ", round(smpd3_ctrl, 1))
  message("  Sgms1 RPKM = ", round(sgms1_ctrl, 1))
  message("  Sgms2 RPKM = ", round(sgms2_ctrl, 1))

  if (smpd3_ctrl > sgms1_ctrl * 1.5) {
    message("  -> Imbalance Smpd3 > Sgms1 en astrocitos (similar al higado)")
  } else {
    message("  -> Smpd3 no supera claramente a Sgms1 en astrocitos")
    message("    (el imbalance Smpd3 > Sms es especifico del higado)")
  }
}

# d) Via S1P en astrocitos
s1pr3_ctrl <- expr_means %>%
  filter(gene == "S1pr3", group == "shControl") %>%
  pull(RPKM_mean)
asah1_ctrl <- expr_means %>%
  filter(gene == "Asah1", group == "shControl") %>%
  pull(RPKM_mean)
sgpl1_ctrl <- expr_means %>%
  filter(gene == "Sgpl1", group == "shControl") %>%
  pull(RPKM_mean)

if (length(s1pr3_ctrl) == 1) {
  message("\nVia S1P en astrocitos EAE (shControl):")
  message("  S1pr3 RPKM = ", round(s1pr3_ctrl, 1))
  message("  Asah1 RPKM = ", round(asah1_ctrl, 1))
  message("  Sgpl1 RPKM = ", round(sgpl1_ctrl, 1))
  if (s1pr3_ctrl > 5) {
    message("  -> S1pr3 esta expresado en astrocitos EAE (consistente con S1P en MS/EAE)")
  }
}

message("\n=== FIN ===")
message("Archivos generados:")
message("  results/24_GSE139531_genes_RPKM_by_sample.csv")
message("  results/24_GSE139531_genes_RPKM_by_group.csv")
message("  results/24_GSE139531_log2FC_vs_shControl.csv")
message("Log: logs/24_GSE139531_expression_log.txt")
