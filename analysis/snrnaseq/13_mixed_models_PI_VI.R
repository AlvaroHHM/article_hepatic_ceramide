# ============================================================
# 13_mixed_models_PI_VI.R
# Regresión lineal mixta + bootstrap para PI/VI
# Testea si MOL3_EAE crece significativamente (Early→Peak→Late)
# y si otros estados EAE resuelven
# ============================================================

# --- Instalar paquetes si faltan ---
for (pkg in c("lme4", "lmerTest", "boot", "ggplot2", "dplyr", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(boot)
  library(ggplot2)
})

dir.create("results/snrnaseq", showWarnings = FALSE, recursive = TRUE)
dir.create("logs/snrnaseq", showWarnings = FALSE, recursive = TRUE)
sink("logs/snrnaseq/13_mixed_models_PI_VI_log.txt", split = TRUE)
on.exit(sink())

cat("=== MODELOS MIXTOS + BOOTSTRAP PARA PI/VI ===\n")
cat("Fecha:", format(Sys.time()), "\n\n")

# ---- 1. Cargar anotación ----
cat("--- 1. Cargando anotación ---\n")
ann_candidates <- c(
  "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz",
  "../data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
)
ann_file <- ann_candidates[file.exists(ann_candidates)][1]
if (is.na(ann_file)) {
  stop("GSE250589 annotation not found. Run 01_oligodendrocytes_GSE250589_annotation.R first.",
       call. = FALSE)
}
annot <- fread(ann_file,
               data.table = FALSE, check.names = FALSE)
rownames(annot) <- annot[[1]]
annot <- annot[, -1]
cat(sprintf("  %d células totales\n", nrow(annot)))

# Definir stage
annot$stage <- ifelse(grepl("^Early", annot$orig.ident), "Early",
               ifelse(grepl("^Peak", annot$orig.ident), "Peak",
               ifelse(grepl("^Late", annot$orig.ident), "Late", NA)))
annot$stage <- factor(annot$stage, levels = c("Early", "Peak", "Late"))

# Filtrar solo EAE
annot_eae <- annot[annot$model == "EAE" & !is.na(annot$stage), ]
cat(sprintf("  %d células EAE\n", nrow(annot_eae)))
cat(sprintf("  %d ratones EAE\n", length(unique(annot_eae$orig.ident))))

cat("\nDistribución de ratones por stage:\n")
print(table(annot_eae$stage, annot_eae$orig.ident))

# ---- 2. Calcular proporción por ratón × estado ----
cat("\n--- 2. Proporción por ratón × estado ---\n")

estados_clave <- c("MOL3_EAE", "MOL5_6_EAE_A", "MOL5_6_EAE_B",
                   "MOL5_6_EAE_C", "MOL1_2_EAE", "MOL2_Cntrl_A")

# Para cada ratón, contar total de células y células por estado
per_mouse <- annot_eae %>%
  group_by(orig.ident, stage) %>%
  summarise(n_total = n(), .groups = "drop")

for (estado in estados_clave) {
  counts <- annot_eae %>%
    filter(cellTypeFalcao2018 == estado) %>%
    group_by(orig.ident) %>%
    summarise(n_state = n(), .groups = "drop")
  per_mouse[[paste0("n_", estado)]] <- 0
  match_idx <- match(counts$orig.ident, per_mouse$orig.ident)
  per_mouse[[paste0("n_", estado)]][match_idx] <- counts$n_state
}

# Proporciones
for (estado in estados_clave) {
  per_mouse[[paste0("prop_", estado)]] <-
    100 * per_mouse[[paste0("n_", estado)]] / per_mouse$n_total
}

cat("\nResumen de proporciones por stage (media ± sd):\n")
for (estado in estados_clave) {
  resumen <- per_mouse %>%
    group_by(stage) %>%
    summarise(media = mean(.data[[paste0("prop_", estado)]]),
              sd = sd(.data[[paste0("prop_", estado)]]),
              n = n(), .groups = "drop")
  cat(sprintf("\n%s:\n", estado))
  print(resumen)
}

write.csv(per_mouse, "results/snrnaseq/13_per_mouse_proportions.csv", row.names = FALSE)

# ---- 3. Modelo mixto por estado ----
cat("\n--- 3. Modelos mixtos (glmer binomial) ---\n")

modelos <- list()
for (estado in estados_clave) {
  n_state <- per_mouse[[paste0("n_", estado)]]
  n_total <- per_mouse$n_total

  # Modelo binomial con random effect por ratón
  tryCatch({
    fit <- glmer(cbind(n_state, n_total - n_state) ~ stage + (1 | orig.ident),
                 data = per_mouse,
                 family = binomial,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    modelos[[estado]] <- fit

    cat(sprintf("\n%s:\n", estado))
    print(summary(fit)$coefficients)

    # Test de tendencia lineal (Early < Peak < Late)
    # Usando la variable numérica de stage
    per_mouse$stage_num <- as.numeric(per_mouse$stage)
    fit_trend <- glmer(cbind(n_state, n_total - n_state) ~ stage_num + (1 | orig.ident),
                       data = per_mouse, family = binomial,
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    coef_trend <- summary(fit_trend)$coefficients["stage_num", ]
    cat(sprintf("  Test de tendencia lineal: coef = %.4f, p = %.4e\n",
                coef_trend["Estimate"], coef_trend["Pr(>|z|)"]))
  }, error = function(e) {
    cat(sprintf("\n%s: ERROR en el modelo - %s\n", estado, e$message))
  })
}

# ---- 4. Bootstrap de PI y VI por ratón ----
cat("\n--- 4. Bootstrap de PI/VI (1000 resamples) ---\n")

calcular_pi_vi <- function(df, estado) {
  # Proporciones por stage (agregando ratones)
  props <- df %>%
    group_by(stage) %>%
    summarise(n_stage = sum(.data[[paste0("n_", estado)]]),
              n_total_stage = sum(n_total),
              prop = 100 * n_stage / n_total_stage, .groups = "drop")
  p_early <- props$prop[props$stage == "Early"]
  p_peak <- props$prop[props$stage == "Peak"]
  p_late <- props$prop[props$stage == "Late"]

  pi <- ifelse(length(p_peak) > 0 && p_peak > 0, p_late / p_peak, NA)
  vi <- ifelse(length(p_early) > 0 && p_early > 0, (p_peak - p_early) / p_early, NA)

  c(pi = pi, vi = vi,
    early = p_early, peak = p_peak, late = p_late)
}

# Bootstrap por ratón (resamplear ratones dentro de cada stage)
set.seed(42)
bootstrap_pi_vi <- function(estado, n_boot = 1000) {
  resultados <- matrix(NA, nrow = n_boot, ncol = 5)
  colnames(resultados) <- c("pi", "vi", "early", "peak", "late")

  ratones_por_stage <- split(per_mouse, per_mouse$stage)

  for (b in 1:n_boot) {
    # Resamplear ratones dentro de cada stage
    boot_df <- do.call(rbind, lapply(ratones_por_stage, function(d) {
      d[sample(nrow(d), replace = TRUE), ]
    }))
    boot_df$orig.ident <- paste0(boot_df$orig.ident, "_b", b)  # evitar duplicados
    res <- calcular_pi_vi(boot_df, estado)
    resultados[b, ] <- res
  }

  data.frame(
    estado = estado,
    pi_mean = mean(resultados[, "pi"], na.rm = TRUE),
    pi_ci_low = quantile(resultados[, "pi"], 0.025, na.rm = TRUE),
    pi_ci_high = quantile(resultados[, "pi"], 0.975, na.rm = TRUE),
    vi_mean = mean(resultados[, "vi"], na.rm = TRUE),
    vi_ci_low = quantile(resultados[, "vi"], 0.025, na.rm = TRUE),
    vi_ci_high = quantile(resultados[, "vi"], 0.975, na.rm = TRUE)
  )
}

cat("\nBootstrap para cada estado (1000 iteraciones):\n")
boot_results <- list()
for (estado in estados_clave) {
  cat(sprintf("  %s...\n", estado))
  boot_results[[estado]] <- bootstrap_pi_vi(estado, n_boot = 1000)
}
boot_df <- do.call(rbind, boot_results)
print(boot_df)
write.csv(boot_df, "results/snrnaseq/13_bootstrap_PI_VI.csv", row.names = FALSE)

# ---- 5. Comparación MOL3_EAE vs otros estados ----
cat("\n--- 5. Comparación MOL3_EAE vs otros estados ---\n")

# Extraer coeficiente de stage para cada estado
coef_table <- data.frame()
for (estado in names(modelos)) {
  fit <- modelos[[estado]]
  coefs <- summary(fit)$coefficients
  # Coeficiente de Peak vs Early (o el primero que exista)
  if ("stagePeak" %in% rownames(coefs)) {
    coef_peak <- coefs["stagePeak", ]
    coef_table <- rbind(coef_table, data.frame(
      estado = estado,
      coef_Peak_vs_Early = coef_peak["Estimate"],
      p_Peak_vs_Early = coef_peak["Pr(>|z|)"],
      stringsAsFactors = FALSE
    ))
  }
}
cat("\nCoeficientes del modelo (Peak vs Early):\n")
print(coef_table)
write.csv(coef_table, "results/snrnaseq/13_coef_table.csv", row.names = FALSE)

# ---- 6. Figura (forest plot de PI/VI con CI) ----
cat("\n--- 6. Generando figura ---\n")

plot_df <- boot_df
plot_df$estado <- factor(plot_df$estado, levels = estados_clave)

# PI plot
p_pi <- ggplot(plot_df, aes(x = estado, y = pi_mean)) +
  geom_point(size = 3, color = "#B2182B") +
  geom_errorbar(aes(ymin = pi_ci_low, ymax = pi_ci_high),
                width = 0.2, color = "#B2182B") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(title = "Persistence Index (PI) con IC 95% (bootstrap, 1000 resamples)",
       subtitle = "PI = %Late / %Peak. Línea discontinua = PI = 1 (sin cambio)",
       x = NULL, y = "PI") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# VI plot
p_vi <- ggplot(plot_df, aes(x = estado, y = vi_mean)) +
  geom_point(size = 3, color = "#2166AC") +
  geom_errorbar(aes(ymin = vi_ci_low, ymax = vi_ci_high),
                width = 0.2, color = "#2166AC") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Velocity Index (VI) con IC 95% (bootstrap, 1000 resamples)",
       subtitle = "VI = (%Peak - %Early) / %Early. Línea discontinua = VI = 0 (sin cambio)",
       x = NULL, y = "VI") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- p_pi / p_vi
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  combined <- arrangeGrob(p_pi, p_vi, nrow = 2)
} else {
  combined <- p_pi
}

ggsave("results/snrnaseq/13_PI_VI_bootstrap.pdf", combined, width = 8, height = 10)
ggsave("results/snrnaseq/13_PI_VI_bootstrap.png", combined, width = 8, height = 10, dpi = 200)
cat("[OK] results/snrnaseq/13_PI_VI_bootstrap.pdf\n")
dir.create("figures/supplementary", showWarnings = FALSE)
file.copy("results/snrnaseq/13_PI_VI_bootstrap.pdf",
          "figures/supplementary/S_bootstrap_PI_VI.pdf", overwrite = TRUE)
cat("[OK] figures/supplementary/S_bootstrap_PI_VI.pdf\n")

# ---- 7. Resumen final ----
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("RESUMEN FINAL\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n", sep = "")

cat("PI y VI con IC 95% (bootstrap):\n")
print(boot_df, row.names = FALSE)

cat("\nCoeficientes de los modelos mixtos:\n")
print(coef_table, row.names = FALSE)

cat("\nINTERPRETACIÓN:\n")
cat("  - Si el IC 95% de PI NO incluye 1 → el cambio es significativo\n")
cat("  - Si el IC 95% de VI NO incluye 0 → el crecimiento es significativo\n")
cat("  - Los estados con IC ancho son los de denominador pequeño\n")
cat("  - MOL3_EAE debe tener IC estrecho si el hallazgo es robusto\n")

cat("\n[OK] Script completado\n")
