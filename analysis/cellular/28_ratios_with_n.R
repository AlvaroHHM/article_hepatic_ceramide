suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
setwd("~/HARDVARD/recount3_analysis")

f <- "results/27_cellular_comparison_long.csv"
if (!file.exists(f)) stop("Missing: ", f)

d <- read.csv(f)
cat("[OK] Filas:", nrow(d), "\n")
cat("[OK] Valores 'cell':", paste(unique(d$cell), collapse = ", "), "\n")

# Pivot SMPD3/SGMS1
sub <- d %>% filter(gene %in% c("Smpd3", "Sgms1"))
piv <- sub %>%
  select(gene, cell, mean_Healthy) %>%
  pivot_wider(names_from = gene, values_from = mean_Healthy) %>%
  mutate(ratio_SMPD3_SGMS1 = Smpd3 / Sgms1)

# n biológico de Methods (n=3 astro, n=8 micro, n=5 neuron, n=5 glia)
n_map <- c(Astrocytes = 3, Microglia = 8, Neurons_POS = 5, Glia_NEG = 5)
piv$n_biological <- n_map[piv$cell]

# Nota ética: descriptive
piv$note <- ifelse(piv$n_biological <= 3,
  "n=3: descriptive; CI not computed",
  "n=5-8: descriptive; not a primary outcome")

cat("\n=== SMPD3 / SGMS1 ratios ===\n")
print(piv[, c("cell", "Smpd3", "Sgms1", "ratio_SMPD3_SGMS1", "n_biological")])

# S1PR1 / S1PR3
piv2 <- d %>% filter(gene %in% c("S1pr1", "S1pr3")) %>%
  select(gene, cell, mean_Healthy) %>%
  pivot_wider(names_from = gene, values_from = mean_Healthy) %>%
  mutate(ratio_S1PR1_S1PR3 = S1pr1 / S1pr3)
piv2$n_biological <- n_map[piv2$cell]

cat("\n=== S1PR1 / S1PR3 ratios ===\n")
print(piv2[, c("cell", "S1pr1", "S1pr3", "ratio_S1PR1_S1PR3", "n_biological")])

# Combinar
piv$ratio_type <- "SMPD3/SGMS1"
piv$ratio_value <- piv$ratio_SMPD3_SGMS1
piv2$ratio_type <- "S1PR1/S1PR3"
piv2$ratio_value <- piv2$ratio_S1PR1_S1PR3
out <- rbind(
  piv[, c("cell", "ratio_type", "ratio_value", "n_biological", "note")],
  data.frame(cell = piv2$cell, ratio_type = piv2$ratio_type,
             ratio_value = piv2$ratio_S1PR1_S1PR3,
             n_biological = piv2$n_biological,
             note = "descriptive; not a primary outcome",
             stringsAsFactors = FALSE)
)

dir.create("results/cellular", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "results/cellular/27b_ratios_with_n.csv", row.names = FALSE)
cat("\n[OK] results/cellular/27b_ratios_with_n.csv\n")
