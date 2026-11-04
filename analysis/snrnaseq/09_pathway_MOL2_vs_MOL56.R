# ============================================================
# 09_pathway_MOL2_vs_MOL56.R
# Comparación de rutas completas entre MOL2_like y MOL56_like
# - Módulo de esfingolípidos (síntesis, interconversión, LacCer, S1P)
# - Módulo de mielina (MBP, PLP1, MAG, MOG, CNP)
# - Módulo de apoptosis (referencia)
# Más robusto al confound que genes individuales
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(hdf5r)
  library(ggplot2); library(dplyr)
})

dir.create("results/snrnaseq", showWarnings = FALSE, recursive = TRUE)
dir.create("logs/snrnaseq", showWarnings = FALSE, recursive = TRUE)
sink("logs/snrnaseq/09_pathway_log.txt", split = TRUE)
on.exit(sink())

cat("=== ANÁLISIS DE RUTA: MOL2_like vs MOL56_like ===\n")
cat("Fecha:", format(Sys.time()), "\n\n")

# ---- 1. Definir módulos de genes ----
pathways <- list(
  "Sphingolipid_synthesis" = c("Sptlc1","Sptlc2","Sptlc3","Kdsr","Degs1","Degs2",
                                "Cers1","Cers2","Cers3","Cers4","Cers5","Cers6"),
  "Sphingolipid_interconversion" = c("Smpd1","Smpd2","Smpd3","Smpd4",
                                      "Sgms1","Sgms2","Asah1","Asah2",
                                      "Acer1","Acer2","Acer3","Sphk1","Sphk2","Sgpl1"),
  "LacCer_synthesis" = c("Ugcg","B4galt5","B4galt6"),
  "S1P_receptors" = c("S1pr1","S1pr2","S1pr3","S1pr4","S1pr5"),
  "Myelin_structural" = c("Mbp","Plp1","Mag","Mog","Cnp"),
  "Apoptosis" = c("Bax","Bak1","Casp3","Casp8","Casp9","Bcl2","Bcl2l1","Bid")
)

cat("--- Módulos definidos ---\n")
for (p in names(pathways)) cat(sprintf("  %-30s: %d genes\n", p, length(pathways[[p]])))
cat("\n")

# ---- 2. Cargar anotación ----
ann_candidates <- c(
  "data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz",
  "../data/GSE250589/GSE250589_EAE_multiomics_annotation.csv.gz"
)
ann_file <- ann_candidates[file.exists(ann_candidates)][1]
if (is.na(ann_file)) {
  stop("GSE250589 annotation not found.\n",
       "To fix: Rscript analysis/snrnaseq/01_GSE250589_annotation.R\n",
       "(downloads 139 MB from GEO into data/GSE250589/)")
}
annot <- fread(ann_file, data.table = FALSE, check.names = FALSE)
rownames(annot) <- annot[[1]]; annot <- annot[, -1]

mol3_bcs_all <- rownames(annot)[annot$cellTypeFalcao2018 == "MOL3_EAE"]
mol3_all <- annot[mol3_bcs_all, ]
mol3_all$subpop <- ifelse(mol3_all$cellType == "MOL2", "MOL2_like",
                   ifelse(mol3_all$cellType == "MOL56", "MOL56_like", "Other"))
mol3_all$stage <- ifelse(grepl("^Early", mol3_all$orig.ident), "Early",
                  ifelse(grepl("^Peak", mol3_all$orig.ident), "Peak",
                  ifelse(grepl("^Late", mol3_all$orig.ident), "Late", "Other")))
mol3_all$barcode <- rownames(mol3_all)

# Subsample estratificado (cap 400) como en 48_v2
set.seed(42)
sampled_bcs <- mol3_all %>%
  slice_sample(prop = 1) %>%
  group_by(subpop, stage) %>%
  slice(1:min(400, n())) %>%
  ungroup() %>% pull(barcode)
mol3_annot <- mol3_all[sampled_bcs, ]
cat(sprintf("Subsample: %d células\n", nrow(mol3_annot)))
print(table(mol3_annot$subpop, mol3_annot$stage))

# ---- 3. Cargar counts de H5 ----
cat("\n--- Cargando counts ---\n")
gsm_map_candidates <- c(
  "results/snrnaseq/00_GSM_mapping.csv",
  "results/snrnaseq/34_GSM_mapping_v2.csv",
  "../results/34_GSM_mapping_v2.csv"
)
gsm_file <- gsm_map_candidates[file.exists(gsm_map_candidates)][1]
if (is.na(gsm_file)) {
  stop("GSM mapping not found. Expected results/snrnaseq/00_GSM_mapping.csv.\n",
       "This file maps orig.ident to GSM in GSE250589.\n",
       "See README 'H5 workflow' section for provenance.")
}
gsm_map <- read.csv(gsm_file, stringsAsFactors = FALSE)
mol3_annot$GSM <- gsm_map$GSM[match(mol3_annot$orig.ident, gsm_map$orig.ident)]

h5_dir <- "data/GSE250589/h5"
if (!dir.exists(h5_dir) || length(list.files(h5_dir, pattern = "\\.h5$")) == 0) {
  stop("H5 files not found at ", h5_dir, "/.\n",
       "To fix: Rscript analysis/snrnaseq/00_download_h5.sh\n",
       "(downloads 28 H5 files, ~2.3 GB, from GEO into data/GSE250589/h5/)")
}
genes_all <- unique(unlist(pathways))
gsms <- unique(mol3_annot$GSM[!is.na(mol3_annot$GSM)])
common_genes <- NULL; h5_paths <- list()
for (gsm in gsms) {
  h5_files <- list.files("data/GSE250589/h5", pattern = gsm, full.names = TRUE)
  if (length(h5_files) == 0) next
  h5 <- H5File$new(h5_files[1], mode = "r")
  g <- as.character(h5[["matrix/features/name"]][])
  h5$close_all()
  h5_paths[[gsm]] <- h5_files[1]
  if (is.null(common_genes)) common_genes <- g
  else common_genes <- intersect(common_genes, g)
}

load_subset <- function(h5_path, wanted_bcs, common_genes) {
  h5 <- H5File$new(h5_path, mode = "r")
  on.exit(tryCatch(h5$close_all(), error = function(e) NULL))
  bcs <- as.character(h5[["matrix/barcodes"]][])
  g <- as.character(h5[["matrix/features/name"]][])
  sh <- as.integer(h5[["matrix/shape"]][])
  d <- as.numeric(h5[["matrix/data"]][])
  i <- as.integer(h5[["matrix/indices"]][])
  p <- as.integer(h5[["matrix/indptr"]][])
  m <- sparseMatrix(i = i + 1L, p = p, x = d, dims = sh, index1 = TRUE,
                   giveCsparse = TRUE)
  rownames(m) <- g; colnames(m) <- bcs
  mi <- match(wanted_bcs, bcs)
  if (all(is.na(mi))) mi <- match(sub("-.*$","",wanted_bcs), sub("-.*$","",bcs))
  valid <- !is.na(mi)
  if (sum(valid) == 0) return(NULL)
  sub <- m[, mi[valid], drop = FALSE]
  colnames(sub) <- wanted_bcs[valid]
  missing <- setdiff(common_genes, rownames(sub))
  if (length(missing) > 0) {
    mm <- Matrix(0, nrow = length(missing), ncol = ncol(sub), sparse = TRUE,
                dimnames = list(missing, colnames(sub)))
    sub <- rbind(sub, mm)
  }
  sub[common_genes, , drop = FALSE]
}

all_mats <- list()
for (gsm in names(h5_paths)) {
  wanted <- rownames(mol3_annot)[mol3_annot$GSM == gsm]
  if (length(wanted) == 0) next
  m <- tryCatch(load_subset(h5_paths[[gsm]], wanted, common_genes),
                error = function(e) NULL)
  if (!is.null(m) && ncol(m) > 0) all_mats[[gsm]] <- m
}

full_mat <- do.call(cbind, all_mats)
mol3_annot <- mol3_annot[colnames(full_mat), ]
cat(sprintf("Matriz: %d genes × %d células\n", nrow(full_mat), ncol(full_mat)))

# ---- 4. Calcular module score ----
cat("\n--- Calculando module scores ---\n")
cpm <- t(t(full_mat) / colSums(full_mat) * 1e6)
logcpm <- log1p(cpm)

# Module score = media de logCPM de los genes de la ruta (presentes en la matriz)
module_scores <- data.frame(barcode = colnames(logcpm))
for (p in names(pathways)) {
  genes_present <- intersect(pathways[[p]], rownames(logcpm))
  if (length(genes_present) < 2) {
    cat(sprintf("  %-30s: SKIP (<2 genes)\n", p))
    next
  }
  module_scores[[p]] <- colMeans(logcpm[genes_present, , drop = FALSE])
  cat(sprintf("  %-30s: %d/%d genes\n", p, length(genes_present), length(pathways[[p]])))
}

# También incluir B4GALT6 solo para comparación
if ("B4galt6" %in% rownames(logcpm)) {
  module_scores[["B4GALT6_alone"]] <- as.numeric(logcpm["B4galt6", ])
}

module_scores$subpop <- mol3_annot$subpop
module_scores$stage <- mol3_annot$stage
module_scores$nCount <- mol3_annot$nCount_RNA

# ---- 5. Análisis estadístico ----
cat("\n--- Comparación por módulo ---\n")

results <- data.frame()
for (p in setdiff(colnames(module_scores), c("barcode","subpop","stage","nCount"))) {
  # Solo MOL2 vs MOL56
  sub <- module_scores[module_scores$subpop %in% c("MOL2_like","MOL56_like"), ]
  w <- wilcox.test(sub[[p]] ~ sub$subpop)
  mol2_mean <- mean(sub[[p]][sub$subpop == "MOL2_like"])
  mol56_mean <- mean(sub[[p]][sub$subpop == "MOL56_like"])
  delta <- mol2_mean - mol56_mean
  pct_delta <- 100 * delta / abs(mol56_mean)
  results <- rbind(results, data.frame(
    module = p,
    MOL2_mean = mol2_mean,
    MOL56_mean = mol56_mean,
    delta = delta,
    pct_delta = pct_delta,
    wilcoxon_p = w$p.value
  ))
}

results <- results[order(-abs(results$delta)), ]
cat("\n")
print(results, row.names = FALSE)
write.csv(results, "results/snrnaseq/09_pathway_MOL2_vs_MOL56.csv", row.names = FALSE)

# ---- 6. Regresión para controlar por depth ----
cat("\n--- Regresión controlando por nCount ---\n")
regression_results <- data.frame()
for (p in setdiff(colnames(module_scores), c("barcode","subpop","stage","nCount"))) {
  sub <- module_scores[module_scores$subpop %in% c("MOL2_like","MOL56_like"), ]
  fit <- lm(sub[[p]] ~ sub$subpop + sub$nCount)
  coef_subpop <- summary(fit)$coefficients["sub$subpopMOL56_like", "Estimate"]
  coef_pval  <- summary(fit)$coefficients["sub$subpopMOL56_like", "Pr(>|t|)"]
  regression_results <- rbind(regression_results, data.frame(
    module = p,
    coef_MOL56_vs_MOL2 = coef_subpop,
    pval = coef_pval
  ))
}
cat("\n")
print(regression_results, row.names = FALSE)
write.csv(regression_results, "results/snrnaseq/09_pathway_regression.csv", row.names = FALSE)

# ---- 7. Plots ----
cat("\n--- Generando plots ---\n")
plot_modules <- c("Sphingolipid_synthesis", "Sphingolipid_interconversion",
                  "LacCer_synthesis", "Myelin_structural", "B4GALT6_alone")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  plots <- lapply(plot_modules[plot_modules %in% colnames(module_scores)],
                  function(p) {
                    ggplot(module_scores[module_scores$subpop != "Other", ],
                           aes(x = subpop, y = .data[[p]], fill = subpop)) +
                      geom_violin(alpha = 0.7) +
                      geom_boxplot(width = 0.15, outlier.shape = NA) +
                      scale_fill_manual(values = c("MOL2_like" = "#B2182B",
                                                    "MOL56_like" = "#2166AC")) +
                      labs(title = p, x = NULL, y = "Module score (logCPM mean)") +
                      theme_minimal(base_size = 9) +
                      theme(legend.position = "none",
                            axis.text.x = element_text(angle = 30, hjust = 1))
                  })
  combined <- wrap_plots(plots, ncol = 2)
  ggsave("results/snrnaseq/09_pathway_MOL2_vs_MOL56.pdf", combined, width = 12, height = 8)
ggsave("results/snrnaseq/09_pathway_MOL2_vs_MOL56.png", combined, width = 12, height = 8, dpi = 150)
dir.create("figures/main", showWarnings = FALSE, recursive = TRUE)
file.copy("results/snrnaseq/09_pathway_MOL2_vs_MOL56.pdf",
          "figures/main/F7_pathway_modules.pdf", overwrite = TRUE)
cat("[OK] figures/main/F7_pathway_modules.pdf\n")
}

# ---- 8. Resumen ----
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("RESUMEN\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n", sep = "")

cat("Efecto observado (sin controlar por depth):\n")
print(results[, c("module","MOL2_mean","MOL56_mean","delta","pct_delta","wilcoxon_p")],
      row.names = FALSE)

cat("\nEfecto después de controlar por nCount (regresión):\n")
print(regression_results, row.names = FALSE)

cat("\nINTERPRETACIÓN:\n")
cat("  - Si |coef| después de regresión < |delta| antes → parte del efecto era depth\n")
cat("  - Si |coef| después de regresión ≈ |delta| antes → el efecto es robusto al depth\n")

cat("\n[OK] results/snrnaseq/09_pathway_MOL2_vs_MOL56.csv\n")
cat("[OK] results/snrnaseq/09_pathway_regression.csv\n")
cat("[OK] results/snrnaseq/09_pathway_MOL2_vs_MOL56.pdf\n")
