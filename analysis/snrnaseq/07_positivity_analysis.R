if (dir.exists("article_hepatic_ceramide")) {
  source("article_hepatic_ceramide/analysis/figures/00_figure_style.R")
} else if (file.exists("analysis/figures/00_figure_style.R")) {
  source("analysis/figures/00_figure_style.R")
}
# 36 - positivity + damage score (v3, robusto)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

ckpt_dir <- "results/35b_checkpoint"
if (!dir.exists(ckpt_dir)) stop("No existe ", ckpt_dir, ". Corre 35b v2 primero.")
ckpts <- list.files(ckpt_dir, pattern = "\\.rds$", full.names = TRUE)
message("Checkpoints: ", length(ckpts))

cells_list <- lapply(ckpts, function(c) {
  x <- tryCatch(readRDS(c), error = function(e) NULL)
  if (is.data.frame(x) && nrow(x) > 0) x else NULL
})
cells_list <- cells_list[!sapply(cells_list, is.null)]
if (length(cells_list) == 0) stop("0 checkpoints con datos")
cells_df <- bind_rows(cells_list)
message("Total cells: ", nrow(cells_df))
message("Columnas: ", paste(colnames(cells_df), collapse = ", "))

genes_int <- c("B4galt6", "S1pr1", "Irf7", "B4galt5", "Ugcg",
               "Smpd3", "Sgms1", "Asah1", "Sgpl1")

# --- Positivity ---
long <- cells_df %>% filter(condition == "EAE") %>%
  pivot_longer(cols = all_of(intersect(genes_int, colnames(cells_df))),
               names_to = "gene", values_to = "CPM") %>%
  mutate(positive = CPM > 0)

positivity <- long %>% group_by(gene, timepoint) %>%
  summarise(n_cells = n(),
            n_positive = sum(positive),
            pct_positive = n_positive / n_cells * 100,
            mean_pos_CPM = ifelse(n_positive > 0, mean(CPM[positive]), 0),
            median_pos_CPM = ifelse(n_positive > 0, median(CPM[positive]), 0),
            .groups = "drop") %>%
  arrange(gene, factor(timepoint, levels = c("Early", "Peak", "Late")))

cat("\n=== POSITIVITY ===\n")
print(positivity, right = FALSE)

ep <- positivity %>% filter(timepoint %in% c("Early", "Peak")) %>%
  select(gene, timepoint, pct_positive, mean_pos_CPM, n_positive) %>%
  pivot_wider(names_from = timepoint,
              values_from = c(pct_positive, mean_pos_CPM, n_positive)) %>%
  mutate(delta_pct = pct_positive_Peak - pct_positive_Early,
         delta_magnitude = log2(mean_pos_CPM_Peak + 1) - log2(mean_pos_CPM_Early + 1),
         mechanism = case_when(
           abs(delta_pct) < 5 & abs(delta_magnitude) < 0.3 ~ "sin_cambio",
           delta_pct > 10 & abs(delta_magnitude) < 0.5   ~ "more_cells",
           abs(delta_pct) < 5 & delta_magnitude > 0.5    ~ "mas_intensidad",
           delta_pct > 10 & delta_magnitude > 0.5        ~ "conversion_+_intensidad",
           delta_pct < -10                               ~ "fewer_cells",
           TRUE                                          ~ "otros")) %>%
  arrange(desc(delta_pct))
cat("\n=== MECANISMO ===\n")
print(ep, right = FALSE)

write.csv(positivity, "results/36_positivity.csv", row.names = FALSE)
write.csv(ep, "results/36_mechanism_early_peak.csv", row.names = FALSE)

# --- Damage score ---
ann_file <- "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
ann <- read.csv(gzfile(ann_file), colClasses = "character", check.names = FALSE)
if (colnames(ann)[1] == "") colnames(ann)[1] <- "barcode"
ann$DamageScore <- as.numeric(ann$DamageScore)
ann$ImmuneScore <- as.numeric(ann$ImmuneScore)
ann$InterferonScore <- as.numeric(ann$InterferonScore)

ann_mol3 <- ann %>% filter(cellTypeFalcao2018 == "MOL3_EAE") %>%
  select(orig.ident, DamageScore, ImmuneScore, InterferonScore)

damage_by_tp <- ann_mol3 %>%
  group_by(orig.ident) %>%
  summarise(mean_Damage = mean(DamageScore, na.rm = TRUE),
            mean_Immune = mean(ImmuneScore, na.rm = TRUE),
            mean_Interferon = mean(InterferonScore, na.rm = TRUE),
            n_mol3 = n(), .groups = "drop") %>%
  left_join(cells_df %>% distinct(orig.ident, timepoint, condition),
            by = "orig.ident") %>%
  filter(condition == "EAE") %>%
  group_by(timepoint) %>%
  summarise(mean_Damage = mean(mean_Damage),
            mean_Immune = mean(mean_Immune),
            mean_Interferon = mean(mean_Interferon),
            n_samples = n(), .groups = "drop")

cat("\n=== DAMAGE/IMMUNE/INTERFERON EN MOL3_EAE (EAE) ===\n")
print(damage_by_tp, right = FALSE)
write.csv(damage_by_tp, "results/36_damage_mol3_by_timepoint.csv", row.names = FALSE)

# --- Plot ---
p <- positivity %>%
  mutate(timepoint = factor(timepoint, levels = c("Early", "Peak", "Late"))) %>%
  filter(gene %in% c("B4galt6", "S1pr1", "Irf7", "Smpd3", "Ugcg")) %>%
  ggplot(aes(x = timepoint, y = pct_positive, group = gene, color = gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "B4galt6" = AAAS$red,
    "S1pr1"   = AAAS$blue,
    "Irf7"    = AAAS$green,
    "Smpd3"   = AAAS$purple,
    "Ugcg"    = AAAS$teal
  ), guide = "none") +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  theme_pub(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 10, face = "bold"),
    panel.grid       = element_blank()
  ) +
  labs(x = NULL, y = "% positive cells (CPM > 0)")
ggsave("results/36_positivity_plot.pdf", p, width = 12, height = 3.5)

message("\nOK: results/36_positivity.csv")
message("OK: results/36_mechanism_early_peak.csv")
message("OK: results/36_damage_mol3_by_timepoint.csv")
message("OK: results/36_positivity_plot.pdf")
