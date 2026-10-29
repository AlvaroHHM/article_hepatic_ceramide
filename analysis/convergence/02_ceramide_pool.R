# ============================================================
# 22_ceramide_pool_and_degradation.R
# (b) Pool completo de ceramidas en PR000633 (NAFLD humano)
# (c) Maquinaria sintesis vs degradacion en higado cronico (SRP074136)
# Uso: Rscript scripts/22_ceramide_pool_and_degradation.R
# ============================================================

suppressPackageStartupMessages({
  library(lipidomeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})


# --- Shared figure style ---
if (dir.exists("article_hepatic_ceramide")) {
  source("article_hepatic_ceramide/analysis/figures/00_figure_style.R")
} else if (file.exists("analysis/figures/00_figure_style.R")) {
  source("analysis/figures/00_figure_style.R")
}
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
sink("logs/22_ceramide_pool_log.txt", split = TRUE)
on.exit(sink())

message("=== Pool de ceramidas + maquinaria de degradacion ===")
message("Fecha: ", Sys.time())

# ============================================================
# PARTE B: Pool completo de ceramidas en PR000633
# ============================================================
message("\n########################################")
message("# PARTE B: Pool de ceramidas en NAFLD")
message("########################################")

data(liverlipidome)

# Todas las especies que empiezan con "Cer(" (ceramidas simples, sin HexCer/DH-Cer)
cer_all <- liverlipidome %>%
  filter(grepl("^Cer\\(", Lipid_Name)) %>%
  filter(!is.na(Lipid_Level))

especies <- unique(cer_all$Lipid_Name)
message("Especies Cer(...) encontradas: ", length(especies))
print(especies)

# Estadisticas por especie y grupo
stats_por_especie <- cer_all %>%
  group_by(Lipid_Name, Diagnosis) %>%
  summarise(
    n      = n(),
    mean   = mean(Lipid_Level, na.rm = TRUE),
    median = median(Lipid_Level, na.rm = TRUE),
    .groups = "drop"
  )

# log2FC vs Normal + t-test por especie
normal_by_especie <- cer_all %>%
  filter(Diagnosis == "Normal") %>%
  group_by(Lipid_Name) %>%
  summarise(normal_mean = mean(Lipid_Level, na.rm = TRUE),
            normal_vals = list(Lipid_Level), .groups = "drop")

tt_by_especie <- do.call(rbind, lapply(especies, function(sp) {
  df_sp <- cer_all %>% filter(Lipid_Name == sp)
  normal_v <- df_sp$Lipid_Level[df_sp$Diagnosis == "Normal"]
  if (length(normal_v) < 2) return(NULL)

  res <- data.frame(Lipid_Name = sp, stringsAsFactors = FALSE)
  for (grp in c("Steatosis","NASH","Cirrhosis")) {
    vals <- df_sp$Lipid_Level[df_sp$Diagnosis == grp]
    if (length(vals) >= 2) {
      tt <- tryCatch(t.test(vals, normal_v), error = function(e) NULL)
      res[[paste0("log2FC_", grp)]] <- log2(mean(vals) + 1) - log2(mean(normal_v) + 1)
      res[[paste0("pval_", grp)]]   <- if (is.null(tt)) NA else tt$p.value
    } else {
      res[[paste0("log2FC_", grp)]] <- NA
      res[[paste0("pval_", grp)]]   <- NA
    }
  }
  res
}))

# Ajuste multiple dentro de cada comparacion
for (grp in c("Steatosis","NASH","Cirrhosis")) {
  tt_by_especie[[paste0("padj_", grp)]] <-
    p.adjust(tt_by_especie[[paste0("pval_", grp)]], method = "BH")
}

# Ordenar por log2FC en NASH (la condicion mas relevante)
tt_by_especie <- tt_by_especie[order(tt_by_especie$log2FC_NASH, decreasing = TRUE), ]

message("\n--- Top 15 especies con mayor aumento en NASH ---")
print(head(tt_by_especie[, c("Lipid_Name",
                              "log2FC_Steatosis", "log2FC_NASH", "log2FC_Cirrhosis",
                              "padj_NASH")], 15),
      row.names = FALSE, digits = 3)

message("\n--- Top 15 especies con mayor disminucion en NASH ---")
print(tail(tt_by_especie[, c("Lipid_Name",
                              "log2FC_Steatosis", "log2FC_NASH", "log2FC_Cirrhosis",
                              "padj_NASH")], 15),
      row.names = FALSE, digits = 3)

# Cuantas especies cambian significativamente?
message("\n--- Especies significativas (padj < 0.05) por comparacion ---")
for (grp in c("Steatosis","NASH","Cirrhosis")) {
  padj_col <- paste0("padj_", grp)
  n_sig <- sum(tt_by_especie[[padj_col]] < 0.05, na.rm = TRUE)
  message("  ", grp, ": ", n_sig, " de ", nrow(tt_by_especie), " especies")
}

# Grafico: las top 6 especies con mayor |log2FC| en NASH
top_up <- head(tt_by_especie[order(tt_by_especie$log2FC_NASH, decreasing = TRUE), ], 3)
top_dn <- head(tt_by_especie[order(tt_by_especie$log2FC_NASH), ], 3)
top6 <- rbind(top_up, top_dn)

cer_top <- cer_all %>% filter(Lipid_Name %in% top6$Lipid_Name)

p_pool <- ggplot(cer_top, aes(x = Diagnosis, y = Lipid_Level, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~ Lipid_Name, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Normal" = "#3B4992",
                                "Steatosis" = "#66CCEE",
                                "NASH" = "#EE0000",
                                "Cirrhosis" = "#631879")) +
  theme_pub(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Lipid_Level")
ggsave("results/22_ceramide_pool_top6.pdf", p_pool, width = 11, height = 6)
message("\nOK results/22_ceramide_pool_top6.pdf")

# ============================================================
# PARTE C: Maquinaria sintesis vs degradacion en higado cronico
# ============================================================
message("\n########################################")
message("# PARTE C: Sintesis vs degradacion en higado cronico")
message("########################################")

liver_file <- "results/18_cers_liver_SRP074136.csv"
if (!file.exists(liver_file)) {
  message("ADVERTENCIA: no encontrado ", liver_file)
  message("Ejecuta primero scripts/18_cers_isoforms_liver_vs_snc.R")
} else {
  liver <- read.csv(liver_file, stringsAsFactors = FALSE)

  # Clasificar genes por funcion
  sintesis <- c("Sptlc1","Sptlc2","Sptlc3","Kdsr","Cers1","Cers2","Cers3",
                "Cers4","Cers5","Cers6","Degs1","Degs2")
  degradacion <- c("Smpd1","Smpd2","Smpd3","Smpd4","Asah1","Asah2",
                   "Acer1","Acer2","Acer3")
  laccer <- c("Ugcg","B4galt5","B4galt6")

  liver$funcion <- NA_character_
  liver$funcion[liver$gene %in% sintesis]    <- "Synthesis"
  liver$funcion[liver$gene %in% degradacion] <- "Degradation"
  liver$funcion[liver$gene %in% laccer]      <- "To_LacCer"

  lfc_col <- "log2FC_inflam_vs_ctrl"

  message("\n--- Maquinaria de SINTESIS en higado cronico ---")
  sint_df <- liver %>% filter(funcion == "Synthesis") %>%
    arrange(desc(.data[[lfc_col]])) %>%
    select(gene, all_of(lfc_col), mean_Mdr2.WT, mean_Mdr2.KO)
  print(sint_df, row.names = FALSE, digits = 2)

  message("\n--- Maquinaria de DEGRADACION en higado cronico ---")
  deg_df <- liver %>% filter(funcion == "Degradation") %>%
    arrange(desc(.data[[lfc_col]])) %>%
    select(gene, all_of(lfc_col), mean_Mdr2.WT, mean_Mdr2.KO)
  print(deg_df, row.names = FALSE, digits = 2)

  message("\n--- Maquinaria de CONVERSION a LacCer ---")
  lac_df <- liver %>% filter(funcion == "To_LacCer") %>%
    arrange(desc(.data[[lfc_col]])) %>%
    select(gene, all_of(lfc_col), mean_Mdr2.WT, mean_Mdr2.KO)
  print(lac_df, row.names = FALSE, digits = 2)

  # Resumen por funcion
  message("\n--- Resumen por funcion (media de log2FC) ---")
  resumen_funcion <- liver %>%
    filter(!is.na(funcion)) %>%
    group_by(funcion) %>%
    summarise(
      n_genes        = n(),
      mean_log2FC    = mean(.data[[lfc_col]], na.rm = TRUE),
      median_log2FC  = median(.data[[lfc_col]], na.rm = TRUE),
      max_log2FC     = max(.data[[lfc_col]], na.rm = TRUE),
      gene_max       = gene[which.max(.data[[lfc_col]])],
      .groups = "drop"
    ) %>%
    arrange(desc(mean_log2FC))
  print(resumen_funcion, row.names = FALSE, digits = 2)

  # Grafico: sintesis vs degradacion vs LacCer
  lfc_col_sym <- rlang::sym(lfc_col)
  p_func <- liver %>%
    filter(!is.na(funcion)) %>%
    ggplot(aes(x = reorder(gene, !!lfc_col_sym),
               y = !!lfc_col_sym, fill = funcion)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    facet_grid(funcion ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("Synthesis" = "#3B4992",
                                  "Degradation" = "#EE0000",
                                  "To_LacCer" = "#631879")) +
    theme_pub(base_size = 9) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "log2FC (Mdr2-KO vs WT)")
  ggsave("results/22_liver_synthesis_vs_degradation.pdf", p_func,
         width = 8, height = 9)
  message("\nOK results/22_liver_synthesis_vs_degradation.pdf")

  # Interpretacion
  message("\n--- INTERPRETACION ---")
  mean_sint <- resumen_funcion$mean_log2FC[resumen_funcion$funcion == "Synthesis"]
  mean_deg  <- resumen_funcion$mean_log2FC[resumen_funcion$funcion == "Degradation"]

  message("Media log2FC sintesis:   ", round(mean_sint, 2))
  message("Media log2FC degradacion: ", round(mean_deg, 2))

  if (!is.na(mean_deg) && !is.na(mean_sint)) {
    if (mean_deg > mean_sint + 0.3) {
      message("\n-> DEGRADACION sube mas que SINTESIS.")
      message("   La ceramida hepatica podria venir de degradacion de esfingomielina")
      message("   (via Smpd3) mas que de sintesis de novo (Cers2).")
      message("   Esto es consistente con un modelo donde la ceramida")
      message("   se acumula sin upregulacion de Cers2.")
    } else if (mean_sint > mean_deg + 0.3) {
      message("\n-> SINTESIS sube mas que DEGRADACION.")
      message("   La ceramida vendria de sintesis de novo, no de degradacion.")
    } else {
      message("\n-> Ambas vias cambian de forma similar.")
      message("   No hay evidencia clara de que una domine sobre la otra.")
    }
  }

  write.csv(sint_df, "results/22_liver_synthesis_genes.csv", row.names = FALSE)
  write.csv(deg_df,  "results/22_liver_degradation_genes.csv", row.names = FALSE)
  write.csv(resumen_funcion, "results/22_liver_summary_by_function.csv",
            row.names = FALSE)
}

# ============================================================
# Guardar
# ============================================================
write.csv(tt_by_especie, "results/22_ceramide_pool_all_species.csv",
          row.names = FALSE)
write.csv(stats_por_especie, "results/22_ceramide_pool_stats.csv",
          row.names = FALSE)

message("\n=== Guardado ===")
message("OK results/22_ceramide_pool_all_species.csv")
message("OK results/22_ceramide_pool_stats.csv")
message("OK results/22_ceramide_pool_top6.pdf")

message("\n=== FIN ===")
