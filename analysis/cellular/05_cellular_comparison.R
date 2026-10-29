# ============================================================
# 27_cellular_comparison.R (v2 corregido)
# Consolida astrocitos, microglia, neuronas en tabla de ratios
# Imprime todo a terminal (y guarda en log)
# Uso: Rscript scripts/27_cellular_comparison.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/27_cellular_comparison_log.txt", split = TRUE)
on.exit(sink())

SEP <- "============================================================"
cat(SEP, "\n")
cat("27_cellular_comparison.R\n")
cat("Fecha:", format(Sys.time()), "\n")
cat(SEP, "\n\n")

# --- 1. Astrocitos (valores del log del script 24, hardcoded) ---
astro <- tribble(
  ~gene,    ~mean_Healthy, ~mean_EAE, ~log2FC,
  "B4galt6",   3.1,  NA,  NA,
  "B4galt5",   2.3,  NA,  NA,
  "Ugcg",      8.6,  NA,  NA,
  "Smpd3",     2.2,  NA,  NA,
  "Sgms1",    21.6,  NA,  NA,
  "Sgms2",     1.6,  NA,  NA,
  "Asah1",    34.7,  NA,  NA,
  "Sgpl1",    21.9,  NA,  NA,
  "S1pr1",   117.5,  NA,  NA,
  "S1pr3",     5.8,  NA,  NA,
  "Pla2g4a",  21.6,  NA,  NA
) %>% mutate(cell = "Astrocytes", dataset = "GSE139531")

# --- 2. Microglia (rename mean_Control -> mean_Healthy) ---
micro_path <- "results/25_GSE253318_microglia_genes.csv"
if (!file.exists(micro_path)) stop("Missing: ", micro_path)
micro <- read.csv(micro_path) %>%
  select(gene, mean_Healthy = mean_Control, mean_EAE, log2FC) %>%
  mutate(cell = "Microglia", dataset = "GSE253318")

# --- 3. Neuronas POS ---
neuron_path <- "results/26_GSE279707_pos.csv"
if (!file.exists(neuron_path)) stop("Missing: ", neuron_path)
neuron <- read.csv(neuron_path) %>%
  select(gene, mean_Healthy, mean_EAE, log2FC) %>%
  mutate(cell = "Neurons_POS", dataset = "GSE279707_POS")

# --- 4. Glia NEG ---
glia_path <- "results/26_GSE279707_neg.csv"
if (!file.exists(glia_path)) stop("Missing: ", glia_path)
glia <- read.csv(glia_path) %>%
  select(gene, mean_Healthy, mean_EAE, log2FC) %>%
  mutate(cell = "Glia_NEG", dataset = "GSE279707_NEG")

# --- 5. Consolidar ---
all <- bind_rows(astro, micro, neuron, glia)

# --- 6. Tabla larga ---
summary_tbl <- all %>%
  select(gene, cell, dataset, mean_Healthy, mean_EAE, log2FC) %>%
  arrange(gene, cell)
write.csv(summary_tbl, "results/27_cellular_comparison_long.csv", row.names = FALSE)

# --- 7. Tabla pivote basal ---
pivot_basal <- summary_tbl %>%
  filter(!is.na(mean_Healthy)) %>%
  select(gene, cell, mean_Healthy) %>%
  pivot_wider(names_from = cell, values_from = mean_Healthy)
write.csv(pivot_basal, "results/27_cellular_comparison_basal.csv", row.names = FALSE)

# --- 8. Ratios (corregido: usa mean_Healthy) ---
ratios <- summary_tbl %>%
  filter(!is.na(mean_Healthy)) %>%
  group_by(cell, dataset) %>%
  summarise(
    SMPD3_SGMS1 = mean_Healthy[gene == "Smpd3"] / mean_Healthy[gene == "Sgms1"],
    S1PR1_S1PR3 = mean_Healthy[gene == "S1pr1"] / mean_Healthy[gene == "S1pr3"],
    B4GALT6     = mean_Healthy[gene == "B4galt6"],
    Smpd3       = mean_Healthy[gene == "Smpd3"],
    Sgms1       = mean_Healthy[gene == "Sgms1"],
    S1pr1       = mean_Healthy[gene == "S1pr1"],
    S1pr3       = mean_Healthy[gene == "S1pr3"],
    .groups = "drop"
  ) %>%
  arrange(cell)
write.csv(ratios, "results/27_cellular_ratios.csv", row.names = FALSE)

# --- 9. IMPRIMIR TODO A TERMINAL ---
cat("\n", SEP, "\n")
cat("TABLA BASAL (mean_Healthy, RPKM/FPKM)\n")
cat(SEP, "\n")
print(as.data.frame(pivot_basal))

cat("\n", SEP, "\n")
cat("RATIOS DIAGNOSTICOS (corregidos)\n")
cat(SEP, "\n")
print(as.data.frame(ratios[, c("cell","Smpd3","Sgms1","SMPD3_SGMS1",
                                "S1pr1","S1pr3","S1PR1_S1PR3","B4GALT6")]))

cat("\n", SEP, "\n")
cat("INTERPRETACION\n")
cat(SEP, "\n")

cat("\n1. Ratio SMPD3/SGMS1 (sintesis ceramida vs esfingomielina):\n")
for (i in seq_len(nrow(ratios))) {
  r <- ratios[i, ]
  fen <- if (r$SMPD3_SGMS1 > 1) "SMPD3-dominante" else "SMS-dominante"
  cat(sprintf("   %-12s: %8.3f  [%s]\n", r$cell, r$SMPD3_SGMS1, fen))
}
cat("\n   Astrocitos son la UNICA celula SMS-dominante (ratio << 1).\n")
cat("   Esto los convierte en receptoras de ceramida (no productoras).\n")

cat("\n2. Ratio S1PR1/S1PR3 (receptor S1P dominante):\n")
for (i in seq_len(nrow(ratios))) {
  r <- ratios[i, ]
  fen <- if (r$S1PR1_S1PR3 > 3) "S1PR1-dominante" else "S1PR3 relevante"
  cat(sprintf("   %-12s: %8.2f  [%s]\n", r$cell, r$S1PR1_S1PR3, fen))
}
cat("\n   S1PR1 es dominante en TODAS las celulas del SNC analizadas.\n")
cat("   Fingolimod (FTY720, aprobado para EM) actua sobre S1PR1.\n")

cat("\n3. B4GALT6 basal:\n")
for (i in seq_len(nrow(ratios))) {
  r <- ratios[i, ]
  cat(sprintf("   %-12s: %8.1f\n", r$cell, r$B4GALT6))
}
cat("\n   Neuronas expresan B4GALT6 mucho mas que astrocitos.\n")
cat("   OJO: plataformas distintas, no comparar magnitudes absolutas.\n")

# --- 10. CAVEAT DE COMPOSICION ---
cat("\n", SEP, "\n")
cat("CAVEAT CRITICO: COMPOSICION EN GSE279707\n")
cat(SEP, "\n")
cat("\nLos log2FC de GSE279707 estan CONFUNDIDOS por cambio de composicion:\n")
cat("  Snap25 cae -0.55 (POS) y -1.62 (NEG)\n")
cat("  Map2   cae -1.33 (POS) y -2.54 (NEG)\n")
cat("  Rbfox3 cae -1.34 (POS) y -1.08 (NEG)\n")
cat("  Gfap   sube +0.74 (POS) y cae -1.22 (NEG)\n")
cat("  Irf7   sube +3.19 (POS) y +3.00 (NEG)\n")
cat("\nInterpretacion: la fraccion NEG en EAE esta dominada por infiltrado\n")
cat("inmune (Irf7 x9, Gfap /2.3), diluyendo marcadores residentes. La\n")
cat("fraccion POS tambien pierde contenido neuronal.\n")
cat("\nCONSECUENCIA: los log2FC de GSE279707 NO reflejan regulacion\n")
cat("intrinseca. Usar SOLO los valores basales (Healthy) como referencia\n")
cat("de expresion por celula; NO usar log2FC de GSE279707.\n")

# --- 11. Heatmap basal ---
genes_ok <- Reduce(intersect, list(astro$gene, micro$gene, neuron$gene, glia$gene))
heat_df <- all %>%
  filter(gene %in% genes_ok, !is.na(mean_Healthy)) %>%
  mutate(log_val = log2(mean_Healthy + 1))

if (nrow(heat_df) > 0) {
  p <- ggplot(heat_df, aes(x = cell, y = gene, fill = log_val)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "white", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = median(heat_df$log_val, na.rm = TRUE),
                         name = "log2(RPKM+1)") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank()) +
    labs(x = NULL, y = NULL,
         title = "Sphingolipid genes across CNS cell types (basal)")
  ggsave("results/27_cellular_heatmap.pdf", p, width = 7, height = 6)
  cat("\nOK results/27_cellular_heatmap.pdf\n")
}

cat("\n", SEP, "\n")
cat("FIN\n")
cat(SEP, "\n")
