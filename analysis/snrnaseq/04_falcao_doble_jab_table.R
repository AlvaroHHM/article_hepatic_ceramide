# ============================================================
# 30b_falcao_doble_jab_table.R  (v3 - fix pct_of_EAE unidades)
# Clasifica los 18 estados Falcao en:
#   homeostatico / EAE-OLG-persistente / EAE-OLG-transitorio /
#   progenitor / microglia_contaminante
# Genera la tabla del doble jab para el anteproyecto v8.
#
# Uso: Rscript scripts/30b_falcao_doble_jab_table.R
# Pre: results/30_falcao_prop_matrix.csv
#      results/30_falcao_persistence.csv
#      results/30_falcao_velocity.csv
#      results/30_falcao_entropy.csv
#
# Cambios v3 vs v2:
#   - FIX: pct_of_EAE calculado como fraccion/fraccion*100
#          (antes: fraccion/porcentaje*100 -> daba fraccion)
#   - Check de dominancia ahora usa >= 50 (porcentaje real)
#   - Mensajes [OK]/[WARN] mas informativos
# ============================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

dir.create("logs", showWarnings = FALSE)
sink("logs/30b_falcao_log.txt", split = TRUE)
on.exit(sink())

message("=== FALCAO DOBLE JAB TABLE (v3) ===")
message("Fecha: ", Sys.time())
message("(fix: pct_of_EAE en porcentaje real, no fraccion)")

# ------------------------------------------------------------
# 1. Cargar CSVs pre-computados (script 30)
# ------------------------------------------------------------
req <- c("results/30_falcao_prop_matrix.csv",
         "results/30_falcao_persistence.csv",
         "results/30_falcao_velocity.csv",
         "results/30_falcao_entropy.csv")
miss <- req[!file.exists(req)]
if (length(miss) > 0) stop("Faltan: ", paste(miss, collapse = ", "),
                           "\nCorre primero scripts/30_falcao_quantification.R")

prop <- read.csv("results/30_falcao_prop_matrix.csv", check.names = FALSE)
pers <- read.csv("results/30_falcao_persistence.csv", check.names = FALSE)
vel  <- read.csv("results/30_falcao_velocity.csv",     check.names = FALSE)
ent  <- read.csv("results/30_falcao_entropy.csv",      check.names = FALSE)

message("\nMatriz de proporciones: ", nrow(prop), " estados x ",
        ncol(prop) - 1, " condiciones")
message("Estados: ", paste(prop$state, collapse = ", "))

# ------------------------------------------------------------
# 2. Clasificar los 18 estados por NOMBRE + metrica PI
# ------------------------------------------------------------
classify_state_v3 <- function(state, pi_late_over_peak) {
  if (grepl("Cntrl", state, ignore.case = TRUE))   return("homeostatico")
  if (grepl("^MiGl", state))                       return("microglia_contaminante")
  if (grepl("^(OPC|COP|NFOL|PLC)", state))         return("progenitor")
  if (grepl("EAE", state)) {
    if (!is.na(pi_late_over_peak) && pi_late_over_peak > 1.0) {
      return("EAE-OLG-persistente")
    }
    return("EAE-OLG-transitorio")
  }
  return("sin_clasificar")
}

pi_lookup <- setNames(pers$PI_Late_over_Peak, pers$state)
vi_lookup <- setNames(vel$VI_Early_to_Peak,  vel$state)

prop$cls <- mapply(classify_state_v3, prop$state, pi_lookup[prop$state])
prop$PI  <- pi_lookup[prop$state]
prop$VI  <- vi_lookup[prop$state]

message("\n--- Clasificacion de los 18 estados ---")
print(prop[, c("state","cls","Early_EAE","Peak_EAE","Late_EAE","PI","VI")],
      row.names = FALSE, digits = 4)

# ------------------------------------------------------------
# 3. Compartimento EAE-OLG strict (excluye homeo/progen/microglia)
# ------------------------------------------------------------
eae_strict <- prop[prop$cls %in% c("EAE-OLG-persistente","EAE-OLG-transitorio"), ]

# Sumas en FRACCION (0-1) - referencia interna
sum_eae_early <- sum(eae_strict$Early_EAE, na.rm = TRUE)
sum_eae_peak  <- sum(eae_strict$Peak_EAE,  na.rm = TRUE)
sum_eae_late  <- sum(eae_strict$Late_EAE,  na.rm = TRUE)

tot_eae <- data.frame(
  time    = c("Early","Peak","Late"),
  pct_eae = c(sum_eae_early, sum_eae_peak, sum_eae_late) * 100
)
message("\n--- Compartimento EAE-OLG strict (% del tejido) ---")
print(tot_eae, row.names = FALSE, digits = 2)

# ------------------------------------------------------------
# 4. TABLA 1: solo estados EAE-OLG (candidatos al doble jab)
# ------------------------------------------------------------
eae_tbl <- eae_strict %>%
  mutate(
    cambio_Early_a_Late = Late_EAE - Early_EAE,
    trayectoria = case_when(
      is.na(VI)                                   ~ "sin datos",
      VI > 0.5  & !is.na(PI) & PI > 1.0           ~ "pico-agudo + persistente",
      VI > 0.5  & !is.na(PI) & PI <= 1.0          ~ "pico-agudo, resuelve",
      VI <= 0.5 & !is.na(PI) & PI > 1.0           ~ "acumulacion monotona",
      TRUE                                        ~ "resolucion"
    )
  ) %>%
  select(state, trayectoria, Early_EAE, Peak_EAE, Late_EAE,
         VI, PI, cambio_Early_a_Late) %>%
  arrange(desc(PI))

message("\n========================================")
message("TABLA 1: Estados EAE-OLG (dobles jab candidatos)")
message("========================================")
print(eae_tbl, row.names = FALSE, digits = 4)

# ------------------------------------------------------------
# 5. TABLA 2: resumen por clase
# ------------------------------------------------------------
summary_by_class <- prop %>%
  group_by(cls) %>%
  summarise(
    n_estados     = n(),
    Early_EAE_pct = sum(Early_EAE, na.rm = TRUE) * 100,
    Peak_EAE_pct  = sum(Peak_EAE,  na.rm = TRUE) * 100,
    Late_EAE_pct  = sum(Late_EAE,  na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(factor(cls, levels = c("EAE-OLG-persistente","EAE-OLG-transitorio",
                                  "homeostatico","progenitor",
                                  "microglia_contaminante","sin_clasificar")))

message("\n========================================")
message("TABLA 2: Resumen por clase de estado")
message("========================================")
print(summary_by_class, row.names = FALSE, digits = 2)

# ------------------------------------------------------------
# 6. TABLA 3: MOL3_EAE como % del compartimento EAE-OLG strict
#    FIX v3: fraccion / fraccion * 100 = porcentaje real
# ------------------------------------------------------------
mol3 <- eae_strict[eae_strict$state == "MOL3_EAE", ]
if (nrow(mol3) != 1) stop("MOL3_EAE no encontrado o duplicado en eae_strict")

tbl3 <- data.frame(
  timepoint     = c("Early","Peak","Late"),
  MOL3_pct      = c(mol3$Early_EAE, mol3$Peak_EAE, mol3$Late_EAE) * 100,
  EAE_total_pct = c(sum_eae_early, sum_eae_peak, sum_eae_late) * 100,
  pct_of_EAE    = c(mol3$Early_EAE / sum_eae_early,
                    mol3$Peak_EAE  / sum_eae_peak,
                    mol3$Late_EAE  / sum_eae_late) * 100
)

message("\n========================================")
message("TABLA 3: MOL3_EAE como % del compartimento EAE-OLG")
message("(fix v3: valores en porcentaje real 0-100)")
message("========================================")
print(tbl3, row.names = FALSE, digits = 2)

# ------------------------------------------------------------
# 7. TABLA 4: desglose de Late EAE (tejido bifurcado)
# ------------------------------------------------------------
late_breakdown <- prop %>%
  select(state, cls, Late_EAE) %>%
  mutate(Late_pct = Late_EAE * 100) %>%
  arrange(desc(Late_pct)) %>%
  filter(Late_pct > 0.5) %>%
  select(state, cls, Late_pct)

message("\n========================================")
message("TABLA 4: Desglose de Late EAE (tejido bifurcado)")
message("========================================")
print(late_breakdown, row.names = FALSE, digits = 2)

late_eae_sum   <- sum_eae_late * 100
late_homeo_sum <- sum(prop$Late_EAE[prop$cls == "homeostatico"], na.rm = TRUE) * 100
late_prog_sum  <- sum(prop$Late_EAE[prop$cls == "progenitor"], na.rm = TRUE) * 100
late_micro_sum <- sum(prop$Late_EAE[prop$cls == "microglia_contaminante"], na.rm = TRUE) * 100

message("\nResumen Late EAE:")
message("  EAE-OLG     = ", round(late_eae_sum,   2), " %")
message("  Homeostatico= ", round(late_homeo_sum, 2), " %")
message("  Progenitor  = ", round(late_prog_sum,  2), " %")
message("  Microglia   = ", round(late_micro_sum, 2), " %")
message("  (Total ~ ", round(late_eae_sum + late_homeo_sum +
                              late_prog_sum + late_micro_sum, 2), " %)")

# ------------------------------------------------------------
# 8. Verificaciones
# ------------------------------------------------------------
message("\n========================================")
message("VERIFICACIONES")
message("========================================")

# 8.1 MOL3_EAE es el UNICO EAE-OLG persistente
persistent_states <- prop$state[prop$cls == "EAE-OLG-persistente"]
if (length(persistent_states) == 1 && persistent_states == "MOL3_EAE") {
  message("[OK] MOL3_EAE es el UNICO estado EAE-OLG persistente")
} else {
  message("[WARN] Estados EAE-OLG persistentes: ",
          paste(persistent_states, collapse = ", "))
}

# 8.2 MOL3_EAE crece monotonicamente
if (mol3$Early_EAE < mol3$Peak_EAE && mol3$Peak_EAE < mol3$Late_EAE) {
  message("[OK] MOL3_EAE crece monotonicamente (Early<Peak<Late)")
} else {
  message("[WARN] MOL3_EAE NO es monotonicamente creciente")
}

# 8.3 MOL3_EAE domina el compartimento EAE en Late (FIX v3)
if (tbl3$pct_of_EAE[3] >= 50) {
  message("[OK] MOL3_EAE domina el compartimento EAE en Late (",
          round(tbl3$pct_of_EAE[3], 1), "%)")
} else {
  message("[WARN] MOL3_EAE no domina Late (",
          round(tbl3$pct_of_EAE[3], 1), "%)")
}

# 8.4 Homeostatico U-shape
late_homeo_frac <- sum(prop$Late_EAE[prop$cls == "homeostatico"], na.rm = TRUE)
peak_homeo_frac <- sum(prop$Peak_EAE[prop$cls == "homeostatico"], na.rm = TRUE)
if (late_homeo_frac > 0.30 && peak_homeo_frac < 0.05) {
  message("[OK] Homeostatico muestra U-shape (colapso en Peak, recuperacion en Late)")
} else {
  message("[WARN] Homeostatico no muestra U-shape clasico")
}

# 8.5 Suma de proporciones ~1 por tiempo
for (t in c("Early_EAE","Peak_EAE","Late_EAE")) {
  s <- sum(prop[[t]], na.rm = TRUE)
  if (abs(s - 1) < 0.02) {
    message("[OK] Suma ", t, " = ", round(s, 4), " (~1.0)")
  } else {
    message("[WARN] Suma ", t, " = ", round(s, 4), " (esperado ~1.0)")
  }
}

# 8.6 Sanity check: pct_of_EAE entre 0 y 100
if (all(tbl3$pct_of_EAE >= 0 & tbl3$pct_of_EAE <= 100)) {
  message("[OK] pct_of_EAE en rango [0, 100] (% real)")
} else {
  message("[WARN] pct_of_EAE fuera de rango: ",
          paste(round(tbl3$pct_of_EAE, 2), collapse = ", "))
}

# ------------------------------------------------------------
# 9. Guardar CSVs
# ------------------------------------------------------------
write.csv(prop,             "results/30b_falcao_state_classification.csv", row.names = FALSE)
write.csv(eae_tbl,          "results/30b_falcao_doble_jab_eae_only.csv",   row.names = FALSE)
write.csv(tbl3,             "results/30b_falcao_mol3_share.csv",           row.names = FALSE)
write.csv(summary_by_class, "results/30b_falcao_summary_by_class.csv",     row.names = FALSE)
write.csv(late_breakdown,   "results/30b_falcao_late_breakdown.csv",       row.names = FALSE)
write.csv(ent,              "results/30b_falcao_entropy_clean.csv",        row.names = FALSE)

message("\n=== Guardado ===")
message("OK results/30b_falcao_state_classification.csv")
message("OK results/30b_falcao_doble_jab_eae_only.csv")
message("OK results/30b_falcao_mol3_share.csv")
message("OK results/30b_falcao_summary_by_class.csv")
message("OK results/30b_falcao_late_breakdown.csv")
message("OK results/30b_falcao_entropy_clean.csv")
message("\n=== FIN ===")
