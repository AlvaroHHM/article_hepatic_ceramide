if (dir.exists("article_hepatic_ceramide")) {
  source("article_hepatic_ceramide/analysis/figures/00_figure_style.R")
} else if (file.exists("analysis/figures/00_figure_style.R")) {
  source("analysis/figures/00_figure_style.R")
}
# ============================================================
# 30_falcao_quantification.R (v2)
# Adaptado al formato real: stage, condition, cellTypeFalcao2018, n, prop
# Entrada:  results/29b_GSE250589_OLG_falcao.csv
# Salida:   results/30_falcao_persistence.csv
#           results/30_falcao_velocity.csv
#           results/30_falcao_entropy.csv
#           results/30_falcao_trajectories.pdf
#           results/30_falcao_entropy.pdf
#           logs/30_falcao_log.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

dir.create("logs", showWarnings = FALSE)
sink("logs/30_falcao_log.txt", split = TRUE)
on.exit(sink())

cat("=== FALCAO QUANTIFICATION v2 (GSE250589) ===\n")
cat("Date:", format(Sys.time()), "\n\n")

# --- Cargar ---
falcao_file <- "results/29b_GSE250589_OLG_falcao.csv"
if (!file.exists(falcao_file)) stop("Falta ", falcao_file)
falcao <- read.csv(falcao_file, check.names = FALSE)
cat("Dimensiones:", nrow(falcao), "x", ncol(falcao), "\n")
cat("Columnas:", paste(colnames(falcao), collapse = ", "), "\n\n")

# --- Verificar columnas esperadas ---
need <- c("stage", "condition", "cellTypeFalcao2018", "n", "prop")
miss <- setdiff(need, colnames(falcao))
if (length(miss) > 0) {
  stop("Faltan columnas: ", paste(miss, collapse = ", "),
       "\nEncontradas: ", paste(colnames(falcao), collapse = ", "))
}

# --- Explorar los valores únicos ---
cat("--- Valores únicos ---\n")
cat("stage:     ", paste(unique(falcao$stage), collapse = ", "), "\n")
cat("condition: ", paste(unique(falcao$condition), collapse = ", "), "\n")
cat("cellType:  ", length(unique(falcao$cellTypeFalcao2018)), "estados\n")
print(unique(falcao$cellTypeFalcao2018))
cat("\n")

# --- Verificar si 'prop' es fracción (0-1) o porcentaje (0-100) ---
prop_check <- falcao %>%
  group_by(stage, condition) %>%
  summarise(sum_prop = sum(prop, na.rm = TRUE),
            n_cells  = sum(n, na.rm = TRUE),
            .groups  = "drop")
cat("--- Suma de prop y n por (stage, condition) ---\n")
print(prop_check)
cat("\nSi sum_prop ~ 1 -> fracción; si ~ 100 -> porcentaje.\n\n")

# --- Pivot a formato wide (una fila por estado) ---
falcao_wide <- falcao %>%
  mutate(combo = paste(stage, condition, sep = "_")) %>%
  select(cellTypeFalcao2018, combo, prop, n) %>%
  pivot_wider(
    names_from  = combo,
    values_from = c(prop, n),
    values_fill = 0
  )

cat("--- Estados en formato wide ---\n")
cat("Columnas:", paste(colnames(falcao_wide), collapse = ", "), "\n\n")

# --- Encontrar columnas de prop para Early/Peak/Late x Ctrl/EAE ---
find_prop <- function(pat) {
  hit <- grep(pat, colnames(falcao_wide), value = TRUE, ignore.case = TRUE)
  if (length(hit) == 0) NA_character_ else hit[1]
}
c_early_ctrl <- find_prop("prop_Early_Ctrl|prop_Ctrl_Early|prop_Early_ctrl")
c_early_eae  <- find_prop("prop_Early_EAE|prop_EAE_Early|prop_Early_eae")
c_peak_ctrl  <- find_prop("prop_Peak_Ctrl|prop_Ctrl_Peak|prop_Peak_ctrl")
c_peak_eae   <- find_prop("prop_Peak_EAE|prop_EAE_Peak|prop_Peak_eae")
c_late_eae   <- find_prop("prop_Late_EAE|prop_EAE_Late|prop_Late_eae")
c_late_ctrl  <- find_prop("prop_Late_Ctrl|prop_Ctrl_Late|prop_Late_ctrl")

cat("--- Columnas de prop detectadas ---\n")
for (k in c("Early_Ctrl","Early_EAE","Peak_Ctrl","Peak_EAE",
            "Late_EAE","Late_Ctrl")) {
  val <- get(paste0("c_", tolower(k)))
  cat(sprintf("  %-12s : %s\n", k, ifelse(is.na(val), "NO ENCONTRADA", val)))
}
cat("\n")

# --- Ensamblar matriz de proporciones ---
safe_pull <- function(df, colname) {
  if (is.na(colname)) return(rep(NA_real_, nrow(df)))
  df[[colname]]
}

prop_matrix <- data.frame(
  state      = falcao_wide$cellTypeFalcao2018,
  Early_Ctrl = safe_pull(falcao_wide, c_early_ctrl),
  Early_EAE  = safe_pull(falcao_wide, c_early_eae),
  Peak_Ctrl  = safe_pull(falcao_wide, c_peak_ctrl),
  Peak_EAE   = safe_pull(falcao_wide, c_peak_eae),
  Late_EAE   = safe_pull(falcao_wide, c_late_eae),
  Late_Ctrl  = safe_pull(falcao_wide, c_late_ctrl),
  stringsAsFactors = FALSE
)
cat("--- Matriz de proporciones ---\n")
print(prop_matrix)
cat("\n")

# --- 1. Persistence Index: Late_EAE / Peak_EAE ---
persistence <- prop_matrix %>%
  mutate(PI_Late_over_Peak = ifelse(!is.na(Peak_EAE) & Peak_EAE > 0,
                                    Late_EAE / Peak_EAE, NA_real_),
         PI_class = case_when(
           is.na(PI_Late_over_Peak)  ~ "sin_dato",
           PI_Late_over_Peak > 1.2   ~ "persistent",
           PI_Late_over_Peak < 0.4   ~ "resolving",
           TRUE                      ~ "parcial"
         )) %>%
  arrange(desc(PI_Late_over_Peak))
cat("--- 1. Persistence Index (Late/Peak) ---\n")
print(persistence %>% select(state, Peak_EAE, Late_EAE, PI_Late_over_Peak, PI_class))
cat("\n")

# --- 2. Velocity Index: Early_EAE -> Peak_EAE ---
velocity <- prop_matrix %>%
  mutate(VI_Early_to_Peak = ifelse(!is.na(Early_EAE) & Early_EAE > 0,
                                   (Peak_EAE - Early_EAE) / Early_EAE,
                                   NA_real_),
         VI_class = case_when(
           is.na(VI_Early_to_Peak)   ~ "sin_dato",
           VI_Early_to_Peak > 1.0    ~ "acelera",
           VI_Early_to_Peak > 0      ~ "crece",
           VI_Early_to_Peak > -0.5   ~ "estable",
           TRUE                      ~ "decae"
         )) %>%
  arrange(desc(VI_Early_to_Peak))
cat("--- 2. Velocity Index (Early->Peak) ---\n")
print(velocity %>% select(state, Early_EAE, Peak_EAE, VI_Early_to_Peak, VI_class))
cat("\n")

# --- 3. Entropía de Shannon por condición ---
shannon_H <- function(p) {
  p <- p[!is.na(p) & p > 0]
  if (length(p) == 0) return(0)
  p <- p / sum(p)
  -sum(p * log2(p))
}

entropy_df <- data.frame(
  combo = c("Early_Ctrl","Early_EAE","Peak_Ctrl","Peak_EAE","Late_EAE","Late_Ctrl"),
  H_bits = c(
    shannon_H(prop_matrix$Early_Ctrl),
    shannon_H(prop_matrix$Early_EAE),
    shannon_H(prop_matrix$Peak_Ctrl),
    shannon_H(prop_matrix$Peak_EAE),
    shannon_H(prop_matrix$Late_EAE),
    shannon_H(prop_matrix$Late_Ctrl)
  ),
  n_states = c(
    sum(prop_matrix$Early_Ctrl > 0, na.rm = TRUE),
    sum(prop_matrix$Early_EAE  > 0, na.rm = TRUE),
    sum(prop_matrix$Peak_Ctrl  > 0, na.rm = TRUE),
    sum(prop_matrix$Peak_EAE   > 0, na.rm = TRUE),
    sum(prop_matrix$Late_EAE   > 0, na.rm = TRUE),
    sum(prop_matrix$Late_Ctrl  > 0, na.rm = TRUE)
  )
)
cat("--- 3. Entropía de Shannon ---\n")
print(entropy_df)
cat("\n")

# --- Guardar CSVs ---
write.csv(persistence, "results/30_falcao_persistence.csv", row.names = FALSE)
write.csv(velocity,    "results/30_falcao_velocity.csv",    row.names = FALSE)
write.csv(entropy_df,  "results/30_falcao_entropy.csv",     row.names = FALSE)
write.csv(prop_matrix, "results/30_falcao_prop_matrix.csv", row.names = FALSE)
cat("OK results/30_falcao_*.csv\n\n")

# --- 4. Figura: trajectories por estado ---
long_traj <- prop_matrix %>%
  select(state, Early_EAE, Peak_EAE, Late_EAE) %>%
  pivot_longer(-state, names_to = "timepoint", values_to = "pct") %>%
  filter(!is.na(pct)) %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("Early_EAE","Peak_EAE","Late_EAE"),
                       labels = c("Early EAE","Peak EAE","Late EAE")),
    trajectory = case_when(
      grepl("MOL3", state)         ~ "persistent (MOL3)",
      grepl("Cntrl", state)        ~ "homeostatic",
      grepl("MOL2", state)         ~ "emerging (MOL2)",
      TRUE                         ~ "resolving (MOL5_6_EAE)"
    )
  )

p_traj <- ggplot(long_traj, aes(x = timepoint, y = pct,
                                 group = state, color = trajectory)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2.8, alpha = 0.9) +
  facet_wrap(~ state, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(
    "persistent (MOL3)"        = "#EE0000",
    "resolving (MOL5_6_EAE)"     = "#3B4992",
    "homeostatic"              = "#008B45",
    "emerging (MOL2)"          = "#631879"
  ), name = "Trajectory") +
  theme_pub(base_size = 9) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "% cells per state")

ggsave("results/30_falcao_trajectories.pdf", p_traj, width = 10, height = 6)
cat("OK results/30_falcao_trajectories.pdf\n")

# --- 5. Figura: entropía ---
p_entropy <- ggplot(entropy_df, aes(x = factor(combo,
                                               levels = c("Early_Ctrl","Early_EAE",
                                                          "Peak_Ctrl","Peak_EAE",
                                                          "Late_EAE","Late_Ctrl")),
                                     y = H_bits, fill = combo)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", H_bits)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c(
    Early_Ctrl = "#008B45", Early_EAE = "#66BD63",
    Peak_Ctrl  = "#D9EF8B", Peak_EAE  = "#631879",
    Late_EAE   = "#EE0000", Late_Ctrl = "#74ADD1"
  ), guide = "none") +
  theme_pub(base_size = 10) +
  labs(x = NULL, y = "Entropía de Shannon (bits)")

ggsave("results/30_falcao_entropy.pdf", p_entropy, width = 7.5, height = 4.5)
cat("OK results/30_falcao_entropy.pdf\n")

cat("\n=== Outputs generados ===\n")
cat("  results/30_falcao_persistence.csv\n")
cat("  results/30_falcao_velocity.csv\n")
cat("  results/30_falcao_entropy.csv\n")
cat("  results/30_falcao_prop_matrix.csv\n")
cat("  results/30_falcao_trajectories.pdf\n")
cat("  results/30_falcao_entropy.pdf\n")
cat("\n=== FIN ===\n")
